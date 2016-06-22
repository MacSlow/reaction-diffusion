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
	BufferType& aBuffer;
	BufferType& bBuffer;
	BufferType& aScratch;
	BufferType& bScratch;
	BufferType& feed;
	BufferType& kill;
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
	_a = BufferType ();
	_a.reserve (_width * _height);
	_b = BufferType ();
	_b.reserve (_width * _height);

	_aScratch = BufferType ();
	_aScratch.reserve (_width * _height);
	_bScratch = BufferType ();
	_bScratch.reserve (_width * _height);

	_f = BufferType ();
	_f.reserve (_width * _height);
	_k = BufferType ();
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
				_b[i + j * _width] = .95;
				_bScratch[i + j * _width] = .05;
			}
		}
	}
}

void Buffer::reset ()
{
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

void updateBuffer (unsigned s,
				   unsigned e,
				   unsigned w,
				   BufferType& aBuffer,
				   BufferType& bBuffer,
				   BufferType& aScratch,
				   BufferType& bScratch,
				   BufferType& feed,
				   BufferType& kill,
				   float dA,
				   float dB,
				   float dt)
{
	for (unsigned index = s; index < e; ++index) {
		float a = aBuffer[index];
		float b = bBuffer[index];
		float f = feed[index];
		float k = kill[index];

		float abb = a * b * b;

		//a, b, _dA, _dB, f, k, _dt, laplaceA(), laplaceB()
		//ymm0 256 bits, 8 floats
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
