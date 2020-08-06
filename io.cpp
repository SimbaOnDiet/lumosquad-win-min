//////////////////////////////////////////////////////////////////////////////////
//
// THOR - A Phyciscal Based Lightning Generator
//
// A Program for the (Undergraduate Course) Electromagnetism Presentation
//
// Modified from LumosQuad by 
// Liang Yien, 
// Student,
// University of Science and Technology of China
// 
///////////////////////////////////////////////////////////////////////////////////
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  I makes no representations about the suitability of this software for
//  any purpose. It is provided "as is" without express or implied warranty.
//  
//  Permission to use, copy, modify and distribute this software and its
//  documentation for educational, research and non-profit purposes, without
//  fee, and without a written agreement is hereby granted, provided that the
//  above copyright notice and the following three paragraphs appear in all
//  copies.
//
//  I SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
//  PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN  "AS IS" BASIS, AND 
//  I HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS,
//  OR MODIFICATIONS.
//
///////////////////////////////////////////////////////////////////////////////////
//
//  This program is largely based on LumosQuad, which has the following copyright 
//  notice and three other free software copyright notice which LumosQuad uses.
//  However, for the reason that I would only like to keep the basic function of
//  generating lightning and make it easy to make others understand the source of 
//  this program, some of them are not included in this program. But I will just
//  keep their copyright notice anyway, just to thanks for their contribution to
//  the original LumosQuad.
//
//  I express my sincere  gratitude to Theodore Kim and Ming C. Lin, the author 
//  of LumosQuad,a wonderful software, and Fast Animation of Lightning Using an 
//  Adaptive Mesh, an inspiring paper, without which I could not have complete 
//  this program.
//  
//  This program also use Mersenne Twister Algorithm to generate random numbers,
//  which is attributed to Makoto Matsumoto and Takuji Nishimura. Its copyright 
//  notice is shown below.
//  
//  
///////////////////////////////////////////////////////////////////////////////////
//
// LumosQuad - A Lightning Generator
// Copyright 2007
// The University of North Carolina at Chapel Hill
//
//  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY WARRANTIES,
//  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
//  FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN
//  "AS IS" BASIS, AND THE UNIVERSITY OF NORTH CAROLINA HAS NO OBLIGATION TO
//  PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
//  Please send questions and comments about LumosQuad to kim@cs.unc.edu.
// 
///////////////////////////////////////////////////////////////////////////////////
//  OpenEXR Copyright Notice:
//
//  Copyright (c) 2002, Industrial Light & Magic, a division of Lucas
//  Digital Ltd. LLC
// 
//  All rights reserved.
// 
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//  *       Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  *       Redistributions in binary form must reproduce the above
//  copyright notice, this list of conditions and the following disclaimer
//  in the documentation and/or other materials provided with the
//  distribution.
//  *       Neither the name of Industrial Light & Magic nor the names of
//  its contributors may be used to endorse or promote products derived
//  from this software without specific prior written permission. 
// 
////////////////////////////////////////////////////////////////////////////////////
//
// This program contains a very thin wrapper to Daniel Dunbar's blue noise generator.
//
// For the original, untainted code, see: 
//   http://www.cs.virginia.edu/~gfx/pubs/antimony/
//
///////////////////////////////////////////////////////////////////////////////////
// This program used Mersenne Twister Random Number Generator, witten by 
// Takuji Nishimura and Makoto Matsumoto, and modified to be a C++ class 
// by Daniel Dunbar. Its copyright notice is shown below:
//
//   A C-program for MT19937, with initialization improved 2002/1/26.
//   Coded by Takuji Nishimura and Makoto Matsumoto.
//   Modified to be a C++ class by Daniel Dunbar.
//
//   Before using, initialize the state by using init_genrand(seed)  
//   or init_by_array(init_key, key_length).
//
//  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//  All rights reserved.                          
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     3. The names of its contributors may not be used to endorse or promote 
//        products derived from this software without specific prior written 
//        permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//   Any feedback is very welcome.
//   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
//
//////////////////////////////////////////////////////////////////////////
//  This program also uses FFTW. Here`s its copyright notice:
//
//
// Copyright (c) 2003, 2007-14 Matteo Frigo
// Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// FFTW, online at http://www.fftw.org
//
////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <iostream>
#include "param.h"
#include <string>
#include "quad_poisson_2d.h"
#include "cell.h"
#include "ppm.h"
#include "Quad_Poisson.h"

static QUAD_DBM_2D* potential = new QUAD_DBM_2D(x_Res, y_Res, Iterations);

bool loadImages(std::string inputFile,int inputWidth,int inputHeight)
{
	// load the files
	unsigned char* input = NULL;
	LoadPPM(inputFile.c_str(), input, inputWidth, inputHeight);

	unsigned char* start = new unsigned char[inputWidth * inputHeight];
	unsigned char* repulsor = new unsigned char[inputWidth * inputHeight];
	unsigned char* attractor = new unsigned char[inputWidth * inputHeight];
	unsigned char* terminators = new unsigned char[inputWidth * inputHeight];

	// composite RGB channels into one
	for (int x = 0; x < inputWidth * inputHeight; x++)
	{
		start[x] = (input[3 * x] == 255) ? 255 : 0;
		repulsor[x] = (input[3 * x + 1] == 255) ? 255 : 0;
		attractor[x] = (input[3 * x + 2] == 255) ? 255 : 0;
		terminators[x] = 0;

		if (input[3 * x] + input[3 * x + 1] + input[3 * x + 2] == 255 * 3)
		{
			terminators[x] = 255;
			start[x] = repulsor[x] = attractor[x] = 0;
		}
	}

	if (potential) delete potential;
	potential = new QUAD_DBM_2D(inputWidth, inputHeight, Iterations);
	bool success = potential->readImage(start, attractor, repulsor, terminators, inputWidth, inputHeight);
	if (success) {
		cout << "File accepted!" << endl;
		cout << "Now working on it." << endl;
		while (potential->addParticle()) {
			//GO! Loop until hit ground.
		}
	}
	else {
		cout << "File rejected!" << endl;
	}
    
	// delete the memory
	delete[] input;
	delete[] start;
	delete[] repulsor;
	delete[] attractor;
	delete[] terminators;

	return success;
}

