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

#ifndef _BUFFER_H
#define _BUFFER_H

#include <vector>
#include <tuple>
#include <memory>
#include <SDL.h>

typedef std::vector<float> BufferType;

class Buffer
{
	public:
		Buffer (unsigned width,
				unsigned height,
				float dA,
				float dB,
				float feed,
				float kill);
		~Buffer ();

		void seed (unsigned x, unsigned y, unsigned ri);
		void reset ();
		void update ();
		void updateMT ();
		void paint (SDL_Window* window);

	private:
		unsigned _width;
		unsigned _height;
		float _dt;
		float _dA;
		float _dB;
		float _feed;
		float _kill;
		BufferType _a;
		BufferType _b;
		BufferType _aScratch;
		BufferType _bScratch;
		BufferType _f;
		BufferType _k;
};

#endif // _BUFFER_H
