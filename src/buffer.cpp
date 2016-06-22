////////////////////////////////////////////////////////////////////////////////
//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
//
// reaction-diffusion simulation in 2D
//
// Copyright 2016 Mirco Müller
//
// Author(s):
//   Mirco "MacSlow" Müller <macslow@gmail.com>
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3, as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranties of
// MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <future>

#include "buffer.h"

extern "C" float asmLaplaceSSE (float* buffer,
								unsigned offset,
								unsigned width,
								unsigned size);

using namespace std;

#define NUM_CHANNELS 4
#define MAX_THREADS 4

struct Input
{
	unsigned start; 
	unsigned end;
	unsigned width;
	BufferType& buffer;
	BufferType& scratch;
	BufferType2& aBuffer;
	BufferType2& bBuffer;
	BufferType2& aScratch;
	BufferType2& bScratch;
	BufferType2& feed;
	BufferType2& kill;
	float dA;
	float dB;
	float dt;
};

Buffer::Buffer (unsigned width,
				unsigned height,
				float dA,
				float dB,
				float feed,
				float kill) :

	_width (width),
	_height (height),
	_dt (.75),
	_dA (dA),
	_dB (dB),
	_feed (feed),
	_kill (kill)
{
	_buffer = vector<tuple<float, float, float, float>> ();
	_buffer.reserve (_width * _height);
	_scratch = vector<tuple<float, float, float, float>> ();
	_scratch.reserve (_width * _height);
	for (unsigned x = 0; x < _width ; ++x) {
		for (unsigned y = 0; y < _height; ++y) {
			tuple<float, float, float, float> t (1.,
												 .0,
												 .01 + ((float) y / (float) _height) * (.1 - .01),
												 .045 + ((float) x / (float) _width) * (.07 - .045));
			_buffer.emplace_back (t);
			_scratch.emplace_back (t);
		}
	}

	_a = BufferType2 ();
	_a.reserve (_width * _height);
	_b = BufferType2 ();
	_b.reserve (_width * _height);

	_aScratch = BufferType2 ();
	_aScratch.reserve (_width * _height);
	_bScratch = BufferType2 ();
	_bScratch.reserve (_width * _height);

	_f = BufferType2 ();
	_f.reserve (_width * _height);
	_k = BufferType2 ();
	_k.reserve (_width * _height);
	for (unsigned x = 0; x < _width ; ++x) {
		for (unsigned y = 0; y < _height; ++y) {
			_a.push_back (1.);
			_b.push_back (.0);
			_aScratch.push_back (1.);
			_bScratch.push_back (.0);
			_f.push_back (.01 + ((float) y / (float) _height) * (.1 - .01));
			_k.push_back (.045 + ((float) x / (float) _width) * (.07 - .045));
		}
	}
}

Buffer::~Buffer ()
{
}

void Buffer::seed (unsigned x, unsigned y, unsigned ri)
{
	float p = .0;
	float q = .0;
	float r = (float) ri;
	for (unsigned i = x - ri; i < x + ri; ++i) {
		for (unsigned j = y - ri; j < y + ri; ++j) {
			p = (float) i - (float) x;
			q = (float) j - (float) y;
			if (p*p + q*q < r*r &&
				i < _width &&
				i >= 1 &&
				j < _height &&
				j >= 1) {
				get<1> (_buffer.at (i + j * _width)) = .95;
				get<1> (_scratch.at (i + j * _width)) = .05;
				_b[i + j * _width] = .95;
				_bScratch[i + j * _width] = .05;
			}
		}
	}
}

void Buffer::reset ()
{
	for (auto& t : _buffer) {
		get<0> (t) = .95;
		get<1> (t) = .05;
	}

	for (auto& t : _scratch) {
		get<0> (t) = .95;
		get<1> (t) = .05;
	}

	for (auto& i : _a) {
		i = .95;
	}
	for (auto& i : _b) {
		i = .05;
	}
	for (auto& i : _aScratch) {
		i = .95;
	}
	for (auto& i : _bScratch) {
		i = .05;
	}
}

float Buffer::laplaceA (unsigned x, unsigned y)
{
	float sum = .0;

	sum += get<0> (_buffer[x + y * _width]) * (-1);
	sum += get<0> (_buffer[x - 1 + y * _width]) * .2;
	sum += get<0> (_buffer[x + 1 + y * _width]) * .2;
	sum += get<0> (_buffer[x + (y - 1) * _width]) * .2;
	sum += get<0> (_buffer[x + (y + 1) * _width]) * .2;
	sum += get<0> (_buffer[x - 1 + (y - 1) * _width]) * .05;
	sum += get<0> (_buffer[x - 1 + (y + 1) * _width]) * .05;
	sum += get<0> (_buffer[x + 1 + (y - 1) * _width]) * .05;
	sum += get<0> (_buffer[x + 1 + (y + 1) * _width]) * .05;

	return sum;
}

float Buffer::laplaceB (unsigned x, unsigned y)
{
	float sum = .0;

	sum += get<1> (_buffer[x + y * _width]) * (-1);
	sum += get<1> (_buffer[x - 1 + y * _width]) * .2;
	sum += get<1> (_buffer[x + 1 + y * _width]) * .2;
	sum += get<1> (_buffer[x + (y - 1) * _width]) * .2;
	sum += get<1> (_buffer[x + (y + 1) * _width]) * .2;
	sum += get<1> (_buffer[x - 1 + (y - 1) * _width]) * .05;
	sum += get<1> (_buffer[x - 1 + (y + 1) * _width]) * .05;
	sum += get<1> (_buffer[x + 1 + (y - 1) * _width]) * .05;
	sum += get<1> (_buffer[x + 1 + (y + 1) * _width]) * .05;

	return sum;
}

float laplaceA3 (BufferType& buffer, unsigned i, unsigned w)
{
	float sum = .0;

	sum += get<0> (buffer[i]) * (-1);
	sum += get<0> (buffer[i - 1]) * .2;
	sum += get<0> (buffer[i + 1]) * .2;
	sum += get<0> (buffer[i - w]) * .2;
	sum += get<0> (buffer[i + w]) * .2;
	sum += get<0> (buffer[i - 1 - w]) * .05;
	sum += get<0> (buffer[i - 1 + w]) * .05;
	sum += get<0> (buffer[i + 1 - w]) * .05;
	sum += get<0> (buffer[i + 1 + w]) * .05;

	return sum;
}

float laplaceB3 (BufferType& buffer, unsigned i, unsigned w)
{
	float sum = .0;

	sum += get<1> (buffer[i]) * (-1);
	sum += get<1> (buffer[i - 1]) * .2;
	sum += get<1> (buffer[i + 1]) * .2;
	sum += get<1> (buffer[i - w]) * .2;
	sum += get<1> (buffer[i + w]) * .2;
	sum += get<1> (buffer[i - 1 - w]) * .05;
	sum += get<1> (buffer[i - 1 + w]) * .05;
	sum += get<1> (buffer[i + 1 - w]) * .05;
	sum += get<1> (buffer[i + 1 + w]) * .05;

	return sum;
}

float laplace (float* buf, unsigned i, unsigned w, unsigned size)
{
	float sum = .0;

	sum += buf[i] * (-1);
	sum += buf[i - 1] * .2;
	sum += buf[i + 1] * .2;
	sum += buf[i - w] * .2;
	sum += buf[i + w] * .2;
	sum += buf[i - 1 - w] * .05;
	sum += buf[i - 1 + w] * .05;
	sum += buf[i + 1 - w] * .05;
	sum += buf[i + 1 + w] * .05;

	return sum;
}

void Buffer::update ()
{
	for (unsigned y = 1; y < _height - 1; ++y) {
		for (unsigned x = 1; x < _width - 1; ++x) {
			float a = get<0> (_buffer[x + y * _width]);
			float b = get<1> (_buffer[x + y * _width]);
			float f = get<2> (_buffer[x + y * _width]);
			float k = get<3> (_buffer[x + y * _width]);
			float abb = a * b * b;

			//a, b, _dA, _dB, f, k, _dt, laplaceA(), laplaceB()
			//ymm0 256 bits, 8 floats
			// x + (y * z - u + v * (1 - w) * r;
			// A = A * B + C (fused multiply-add)
			get<0> (_scratch.at (x + y * _width)) = a + (_dA *
														 laplaceA (x, y) -
														 abb +
														 f * (1. - a)) * _dt;
			get<1> (_scratch.at (x + y * _width)) = b + (_dB *
														 laplaceB (x, y) +
														 abb -
														 (f + k) * b) * _dt;
		}
	}
	swap (_buffer, _scratch);
}

void updateBuffer (unsigned s,
				   unsigned e,
				   unsigned w,
				   BufferType& buffer,
				   BufferType& scratch,
				   BufferType2& aBuffer,
				   BufferType2& bBuffer,
				   BufferType2& aScratch,
				   BufferType2& bScratch,
				   BufferType2& feed,
				   BufferType2& kill,
				   float dA,
				   float dB,
				   float dt)
{
	for (unsigned index = s; index < e; ++index) {
		/*float a = get<0> (buffer[index]);
		float b = get<1> (buffer[index]);
		float f = get<2> (buffer[index]);
		float k = get<3> (buffer[index]);*/
		float a = aBuffer[index];
		float b = bBuffer[index];
		float f = feed[index];
		float k = kill[index];

		float abb = a * b * b;

		//a, b, _dA, _dB, f, k, _dt, laplaceA(), laplaceB()
		//ymm0 256 bits, 8 floats
		/*get<0> (scratch.at (index)) = a + (dA *
										   laplaceA3 (buffer, index, w) -
										   abb +
										   f * (1. - a)) * dt;
		get<1> (scratch.at (index)) = b + (dB *
											laplaceB3 (buffer, index, w) +
											abb -
											(f + k) * b) * dt;*/
		aScratch[index] = a + (dA *
							   //asmLaplaceSSE (aBuffer.data(), index, w, 4) -
							   laplace (aBuffer.data(), index, w, 4) -
							   abb +
							   f * (1. - a)) * dt;
		bScratch[index] = b + (dB *
							   //asmLaplaceSSE (bBuffer.data(), index, w, 4) +
							   laplace (bBuffer.data(), index, w, 4) +
							   abb -
							   (f + k) * b) * dt;
	}
}

void updateBufferThreaded (future<Input>& f)
{
    Input input = f.get ();

    updateBuffer (input.start,
				  input.end,
                  input.width,
                  input.buffer,
                  input.scratch,
                  input.aBuffer,
                  input.bBuffer,
                  input.aScratch,
                  input.bScratch,
                  input.feed,
                  input.kill,
                  input.dA,
                  input.dB,
                  input.dt);
}

void Buffer::updateMT ()
{
	promise<Input> p[MAX_THREADS];
	future<Input> fInput[MAX_THREADS];
	future<void> f[MAX_THREADS];
	Input input = {0,
				   0,
				   _width,
				   _buffer,
				   _scratch,
				   _a,
				   _b,
				   _aScratch,
				   _bScratch,
				   _f,
				   _k,
				   _dA,
				   _dB,
				   _dt};

	unsigned size = _width * (_height - 1) - (_width + 1);
	size /= MAX_THREADS;

	for (unsigned i = 0; i < MAX_THREADS; ++i) {
		fInput[i] = p[i].get_future ();
		f[i] = future<void> (async (launch::async,
				   //updateBufferAVXthreaded,
				   updateBufferThreaded,
				   ref (fInput[i])));
		input.start = _width + 1 + i * size;
		input.end   = _width + 1 + (i+1) * size;
		p[i].set_value (input);
	}

	for (unsigned i = 0; i < MAX_THREADS; ++i) {
		f[i].get ();
	}

	//swap (_buffer, _scratch);
	swap (_a, _aScratch);
	swap (_b, _bScratch);
	/*float buf[] = { 1.,  2.,  3.,  4.,  1.,  2.,  3.,  4.,
				    5.,  6.,  7.,  8.,  5.,  6.,  7.,  8.,
				    9., 10., 11.,  12.,  9., 10., 11., 12.,
	                2.,  4.,  6.,  8.,  2.,  4.,  6.,  8.,
				   10., 12., 14., 16., 10., 12., 14., 16.,
				   18., 20., 22., 24., 18., 20., 22., 24.};
	cout << "asmLaplace(): " << asmLaplaceSSE (buf, 36, 8, 4) << endl;
	cout << "asmLaplace(): " << asmLaplaceSSE (buf, 52, 8, 4) << endl;
	cout << "asmLaplace(): " << asmLaplaceSSE (buf, 100, 8, 4) << endl;
	cout << "asmLaplace(): " << asmLaplaceSSE (buf, 132, 8, 4) << endl;
	cout << "asmLaplace(): " << asmLaplaceSSE (buf, 148, 8, 4) << endl << endl;*/
}

void Buffer::paint (SDL_Window* window)
{
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
	unsigned rmask = 0xff000000;
	unsigned gmask = 0x00ff0000;
	unsigned bmask = 0x0000ff00;
	unsigned amask = 0x000000ff;
#else
	unsigned rmask = 0x000000ff;
	unsigned gmask = 0x0000ff00;
	unsigned bmask = 0x00ff0000;
	unsigned amask = 0xff000000;
#endif
	unsigned pitch = _width * NUM_CHANNELS;
	unsigned size = _height * pitch;
	unsigned char* buffer = (unsigned char*) calloc (size, sizeof (unsigned char));
	SDL_Surface* surface = nullptr;

	unsigned indexA = 0;
	unsigned indexB = 0;
	float value = 0;

	for (unsigned y = 0; y < _height; ++y) {
		for (unsigned x = 0; x < _width; ++x) {
			indexA = x * NUM_CHANNELS + y * pitch;
			indexB = x + y * _width;
			//value = get<0> (_buffer[indexB]);
			value = _a[indexB];
			buffer[indexA] = (Uint8) (value * 255.);
			buffer[indexA+1] = (Uint8) (value * 255.);
			buffer[indexA+2] = (Uint8) (value * 255.);
			buffer[indexA+3] = 255;
		}
	}

	surface = SDL_CreateRGBSurfaceFrom (buffer,
										_width,
										_height,
										NUM_CHANNELS * 8,
										pitch,
										rmask,
										gmask,
										bmask,
										amask);
	SDL_Surface* src = NULL;
	src = SDL_ConvertSurfaceFormat (surface, SDL_PIXELFORMAT_RGB888, 0);
	SDL_FreeSurface (surface);
	free (buffer);

	SDL_Surface* dst = SDL_GetWindowSurface (window);
    SDL_BlitSurface (src, NULL, dst, NULL);
    SDL_FreeSurface (src);
    SDL_UpdateWindowSurface (window);
}
