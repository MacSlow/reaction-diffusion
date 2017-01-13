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

#include <iostream>
#include <sstream>

#include "display.h"

using namespace std;

#define WIN_TITLE "Reaction-Diffusion by MacSlow"

void Display::initialize ()
{
	if (_initialized) {
		return;
	}

	// initialize SDL
	int result = 0;
	SDL_ClearError ();
	result = SDL_Init (SDL_INIT_VIDEO);
	if (result != 0) {
		std::cout << "SDL_Init() failed: " << SDL_GetError () << std::endl;
		_initialized = false;
		return;
	}

    _initialized = true;
}

Display::Display (unsigned int width, unsigned int height) :
	_initialized (false),
	_window (NULL),
	_running (false),
	_paused (false),
	_lmb (false),
	_penSize (4),
	_seconds (0),
	_allFrames (0),
	_min (0),
	_avg (0),
	_max (0)
{
	initialize ();

	SDL_ClearError ();
	_window = SDL_CreateWindow (WIN_TITLE,
								SDL_WINDOWPOS_UNDEFINED,
								SDL_WINDOWPOS_UNDEFINED,
								width,
								height,
								0);
	if (!_window) {
		std::cout << "window creation failed: " << SDL_GetError () << std::endl;
		return;
	}
	_buffer = make_unique<Buffer> (width, height, .975, .2, .078, .061);

	// "dump" some initial substance A into B at different spots
	_buffer->seed (150, 150, 50);
	_buffer->seed (350, 150, 20);
	_buffer->seed (250, 250, 10);
	_buffer->seed (150, 350, 40);
	_buffer->seed (350, 350, 30);
}

Display::~Display ()
{
	SDL_DestroyWindow (_window);
	SDL_Quit ();
}

bool Display::run ()
{
	if (!_initialized) {
		return false;
	}

	_running = true;

	while (_running) {
		SDL_Event event;
		while (SDL_PollEvent (&event)) {
			switch (event.type) {
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == SDLK_RETURN) {
						_buffer->reset ();
					}

					if (event.key.keysym.sym == SDLK_ESCAPE) {
						_running = false;
					} else if (event.key.keysym.sym == SDLK_SPACE) {
						_paused = !_paused;
						if (_paused) {
							std::stringstream title;
							title << WIN_TITLE << " - paused";
							std::string str (title.str ());
							SDL_SetWindowTitle (_window, str.c_str ());
						} else {
							SDL_SetWindowTitle (_window, WIN_TITLE);
						}
					}
				break;

				case SDL_MOUSEBUTTONDOWN:
					int y;
					SDL_GetMouseState (nullptr, &y);
					if (event.button.button == SDL_BUTTON_LEFT && y > 0) {
						_lmb = true;
					}
				break;

				case SDL_MOUSEBUTTONUP:
					if (event.button.button == SDL_BUTTON_LEFT) {
						_lmb = false;
					}
				break;

				case SDL_MOUSEWHEEL:
					if (event.wheel.y == 1) {
						_penSize = _penSize < 30 ? _penSize + 1 : _penSize;
					} else if (event.wheel.y == -1) {
						_penSize = _penSize > 3 ? _penSize - 1 : _penSize;
					}
				break;

				case SDL_MOUSEMOTION:
					if (_lmb) {
						int x;
						int y;
						SDL_GetMouseState (&x, &y);
						_buffer->seed (x, y, _penSize);
					}
				break;

				case SDL_QUIT:
					_running = false;
				break;
			}
		}

		update ();
	}

	return true;
}

bool Display::update ()
{
	if (!_initialized) {
		return false;
	}

	if (_paused) {
		return true;
	}

	static unsigned int fps = 0;
	static unsigned int lastTick = 0;
	static unsigned int currentTick = 0;

	// spit out frame-rate and frame-time
	fps++;
	currentTick = SDL_GetTicks ();
	if (currentTick - lastTick > 1000) {
		std::stringstream title;
		title << WIN_TITLE << " - "
			  << _min << "/"
			  << _avg << "/"
			  << _max << " fps";
		std::string str (title.str ());
		SDL_SetWindowTitle (_window, str.c_str ());
		_seconds++;
		lastTick = currentTick;

		if (_max < fps) {
			_max = fps;
		}
		if (_seconds > 0) {
			_allFrames += fps;
			_avg = _allFrames / _seconds;
		}
		if (_min > fps || _min == 0) {
			_min = fps;
		}

		fps = 0;
	}

	_buffer->updateMT ();
	_buffer->paint (_window);

	return true;
}
