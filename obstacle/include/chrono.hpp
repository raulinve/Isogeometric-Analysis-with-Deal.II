//! An example of use of the C++11 chrono library.

/*
This file illustrate the use of the ctime library,

Luca Formaggia 2005. New version using C++11 (2015)

This file is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This file is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
License along with this file; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef _CHRONO_H_INCLUDE_
#define _CHRONO_H_INCLUDE_
#include<chrono>
#include<iosfwd>
#include <iostream>

namespace Timings{

  
  class Chrono
  {
  public:
    using MyClock=std::chrono::high_resolution_clock;
    using MyTimePoint=std::chrono::time_point<MyClock>;
  private:
    MyTimePoint startTime{MyClock::now()};
    MyTimePoint stopTime{MyClock::now()};
  public:
    //! Outputs time from last start and stop
    friend std::ostream & operator <<(std::ostream &out, Chrono const &c) {
	  auto oldf = out.flags();
	  out << "Elapsed Time = " << std::scientific << c.wallTime() << " microsec ( " << c.wallTime()/1000000 << " sec )"
		  << std::endl;
	  out.flags(oldf);
	  return out;
	};
    //! Starts/reset  counting time
    void start() {
	  startTime = stopTime = MyClock::now();
	};
    //! Stops counting time
    void stop(){
  	  stopTime = MyClock::now();
	};
    //!Outputs wall time between last start and stop (in microseconds)
    double wallTime() const {
	  using namespace std::chrono;
	  auto time_span = duration_cast<nanoseconds>(stopTime - startTime);
	  return time_span.count() / 1000.;
	};
    //!Outputs wall time between last start and now! (in microseconds)
    double wallTimeNow() const {
	  using namespace std::chrono;
	  auto time_span = duration_cast<nanoseconds>(MyClock::now() - startTime);
	  return time_span.count() / 1000.;
	};
  };

}


#endif
