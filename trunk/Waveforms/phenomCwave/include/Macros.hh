/*
Copyright (C) 2001, 2005 R. Balasubramanian & S. Babak

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
/******************  CVS info ************************ 
#define CVSTAGHH "$Name:  $"
#define CVSHEADERHH "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/Macros.hh,v 1.1 2007/09/18 21:50:37 stba Exp $" 
*/

#ifndef MACROSHH
#define MACROSHH
#include <iostream> 

#define LISAWPAssert(assertion,mesg) if(!(assertion))\
{std::cout<< " Error in line " << __LINE__ << " of file " \
<< __FILE__ << ":" << mesg << std::endl; exit(-1);}

//#define spr "     "; 

#define LISAWPLOGGER 1
#ifdef LISAWPLOGGER
#define LISAWPLogger(outstream, expression) (outstream) << (expression);
#endif
#endif // end of MAcroshh

