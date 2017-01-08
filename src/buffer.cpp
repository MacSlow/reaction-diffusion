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

extern "C" float laplaceSSE (float* buffer,   // RDI
							 unsigned offset, // RSI
							 unsigned width); // RDX

extern "C" void calcA (unsigned index,     // RDI
					   unsigned width,     // RSI
					   float* aScratch,    // RDX
					   float* aBuffer,     // RCX
					   float* bBuffer,     // R8
					   float* feed,        // R9
					   float dA,           // xmm0
					   float dt);          // xmm1

extern "C" void calcB (unsigned index,     // RDI
					   unsigned width,     // RSI
					   float* bScratch,    // RDX
					   float* aBuffer,     // RCX
					   float* bBuffer,     // R8
					   float* feed,        // R9
					   float* kill,        // stack
					   float dB,           // xmm0
					   float dt);          // xmm1

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
	_dt (.75f),
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
			_a.push_back (1.f);
			_b.push_back (.0f);
			_aScratch.push_back (1.f);
			_bScratch.push_back (.0f);
			_f.push_back (.01f + ((float) y / (float) _height) * (.1f-.01f));
			_k.push_back (.045f + ((float) x / (float) _width) * (.07f-.045f));
		}
	}

}

Buffer::~Buffer ()
{
}

void Buffer::seed (unsigned x, unsigned y, unsigned ri)
{
	float p = .0f;
	float q = .0f;
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
				_b[i + j * _width] = .95f;
				_bScratch[i + j * _width] = .95f;
			}
		}
	}
}

void Buffer::reset ()
{
	for (auto& i : _a) {
		i = .95f;
	}
	for (auto& i : _b) {
		i = .05f;
	}
	for (auto& i : _aScratch) {
		i = .95f;
	}
	for (auto& i : _bScratch) {
		i = .05f;
	}
}

float laplaceCPP (float* buf, unsigned i, unsigned w)
{
	float sum = .0f;

	sum += buf[i] * (-1.f);
	sum += buf[i - 1] * .2f;
	sum += buf[i + 1] * .2f;
	sum += buf[i - w] * .2f;
	sum += buf[i + w] * .2f;
	sum += buf[i - 1 - w] * .05f;
	sum += buf[i - 1 + w] * .05f;
	sum += buf[i + 1 - w] * .05f;
	sum += buf[i + 1 + w] * .05f;

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

		aScratch[index] = a + (dA *
							   laplaceSSE (aBuffer.data(), index, w) -
							   abb +
							   f * (1.f - a)) * dt;

		bScratch[index] = b + (dB *
							   laplaceSSE (bBuffer.data(), index, w) +
							   abb -
							   (f + k) * b) * dt;

		/*calcB (index,
			   w,
			   bScratch.data (),
			   aBuffer.data (),
			   bBuffer.data (),
			   feed.data (),
			   kill.data (),
			   dB,
			   dt);*/
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
	unsigned char* buffer = (unsigned char*) calloc (size,
													 sizeof (unsigned char));
	SDL_Surface* surface = nullptr;

	unsigned indexA = 0;
	unsigned indexB = 0;
	float value = 0;

	for (unsigned y = 0; y < _height; ++y) {
		for (unsigned x = 0; x < _width; ++x) {
			indexA = x * NUM_CHANNELS + y * pitch;
			indexB = x + y * _width;
			value = _a[indexB];
			buffer[indexA] = (Uint8) (value * 255.f);
			buffer[indexA+1] = (Uint8) (value * 255.f);
			buffer[indexA+2] = (Uint8) (value * 255.f);
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
