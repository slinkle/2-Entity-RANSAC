/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/


#include <opengv/absolute_pose/LineNoncentralAbsoluteAdapter.hpp>

opengv::absolute_pose::LineNoncentralAbsoluteAdapter::LineNoncentralAbsoluteAdapter(
	  const bearingVectors_t & bearingVectors,
      const camCorrespondences_t & camCorrespondences,
      const points_t & points,
      const bearingVectors_t &  linestart_bearingVectors,
      const bearingVectors_t &  lineend_bearingVectors,
      const camCorrespondences_t & line_camCorrespondences,
      const points_t & linestart_points,
      const points_t & lineend_points ,
      const translations_t & camOffsets,
      const rotations_t & camRotations ) :
    AbsoluteLineAdapterBase(),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _linestart_bearingVectors(linestart_bearingVectors),
    _line_camCorrespondences(line_camCorrespondences),
    _linestart_points(linestart_points),
    _lineend_bearingVectors(lineend_bearingVectors),
    _lineend_points(lineend_points),
    _camOffsets(camOffsets),
    _camRotations(camRotations)
{}

opengv::absolute_pose::LineNoncentralAbsoluteAdapter::LineNoncentralAbsoluteAdapter(
	  const bearingVectors_t & bearingVectors,
      const camCorrespondences_t & camCorrespondences,
      const points_t & points,
      const bearingVectors_t &  linestart_bearingVectors,
      const bearingVectors_t &  lineend_bearingVectors,
      const camCorrespondences_t & line_camCorrespondences,
      const points_t & linestart_points,
      const points_t & lineend_points ,
      const translations_t & camOffsets,
      const rotations_t & camRotations,
    const rotation_t & R  ) :
    AbsoluteLineAdapterBase(R),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _linestart_bearingVectors(linestart_bearingVectors),
    _line_camCorrespondences(line_camCorrespondences),
    _linestart_points(linestart_points),
    _lineend_bearingVectors(lineend_bearingVectors),
    _lineend_points(lineend_points),
    _camOffsets(camOffsets),
    _camRotations(camRotations)
{}

opengv::absolute_pose::LineNoncentralAbsoluteAdapter::LineNoncentralAbsoluteAdapter(
	  const bearingVectors_t & bearingVectors,
      const camCorrespondences_t & camCorrespondences,
      const points_t & points,
      const bearingVectors_t &  linestart_bearingVectors,
      const bearingVectors_t &  lineend_bearingVectors,
      const camCorrespondences_t & line_camCorrespondences,
      const points_t & linestart_points,
      const points_t & lineend_points ,
      const translations_t & camOffsets,
      const rotations_t & camRotations,
    const translation_t & t,
    const rotation_t & R ) :
    AbsoluteLineAdapterBase(t,R),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _linestart_bearingVectors(linestart_bearingVectors),
    _line_camCorrespondences(line_camCorrespondences),
    _linestart_points(linestart_points),
    _lineend_bearingVectors(lineend_bearingVectors),
    _lineend_points(lineend_points),
    _camOffsets(camOffsets),
    _camRotations(camRotations)
{}

opengv::absolute_pose::LineNoncentralAbsoluteAdapter::~LineNoncentralAbsoluteAdapter()
{}

opengv::bearingVector_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getBearingVector( size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _bearingVectors[index];
}

double
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getWeight( size_t index ) const
{
  return 1.0;
}

opengv::point_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getPoint( size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _points[index];
}

opengv::translation_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getCamOffset( size_t index ) const
{
  assert(_camCorrespondences[index] < _camOffsets.size());
  return _camOffsets[_camCorrespondences[index]];
}

opengv::rotation_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getCamRotation( size_t index ) const
{
  assert(_camCorrespondences[index] < _camRotations.size());
  return _camRotations[_camCorrespondences[index]];
}

opengv::bearingVector_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLinestartBearingVector( size_t index ) const
{
  assert(index < _linestart_bearingVectors.size());
  return _linestart_bearingVectors[index];
}

opengv::point_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLinestartPoint( size_t index ) const
{
  assert(index < _linestart_bearingVectors.size());
  return _linestart_points[index];
}

opengv::translation_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLineCamOffset( size_t index ) const
{
  assert(_line_camCorrespondences[index] < _camOffsets.size());
  return _camOffsets[_line_camCorrespondences[index]];
}

opengv::rotation_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLineCamRotation( size_t index ) const
{
  assert(_line_camCorrespondences[index] < _camRotations.size());
  return _camRotations[_line_camCorrespondences[index]];
}

opengv::bearingVector_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLineendBearingVector( size_t index ) const
{
  assert(index < _linestart_bearingVectors.size());
  return _lineend_bearingVectors[index];
}

opengv::point_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLineendPoint( size_t index ) const
{
  assert(index < _lineend_bearingVectors.size());
  return _lineend_points[index];
}


size_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getNumberCorrespondences() const
{
  return _bearingVectors.size();
}

size_t
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::
    getLineNumberCorrespondences() const
{
  return _linestart_bearingVectors.size();
}

size_t 
opengv::absolute_pose::LineNoncentralAbsoluteAdapter::getCamid(std::size_t index) const
{
  return _camCorrespondences[index];
}

