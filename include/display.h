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

#ifndef _DISPLAY_H
#define _DISPLAY_H

#include <memory>
#include <SDL.h>

#include "buffer.h"

class Display
{
	public:
		Display (unsigned int width, unsigned int height);
		~Display ();

		bool run ();
		bool update ();

	private:
		void initialize ();

	private:
		bool _initialized;
		SDL_Window* _window;
		bool _running;
		bool _paused;
		std::unique_ptr<Buffer> _buffer;
		bool _lmb;
		unsigned _penSize;
};

#endif // _DISPLAY_H
