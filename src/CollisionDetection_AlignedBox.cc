/*
(***********************************************************************)
(*                                                                     *)
(* The COLLISION DETECTION project                                     *)
(*                                                                     *)
(* Copyright (c) 2020-2021, Davide Stocco and Enrico Bertolazzi        *)
(*                                                                     *)
(* The COLLISION DETECTION project and its components are supplied     *)
(* under the terms of the open source BSD 2-Clause License. The        *)
(* contents of the COLLISION DETECTION project and its components may  *)
(* not be copied or disclosed except in accordance with the terms of   *)
(* the BSD 2-Clause License.                                           *)
(*                                                                     *)
(* URL: https://opensource.org/licenses/BSD-2-Clause                   *)
(*                                                                     *)
(*    Davide Stocco                                                    *)
(*    Department of Industrial Engineering                             *)
(*    University of Trento                                             *)
(*    e-mail: davide.stocco@unitn.it                                   *)
(*                                                                     *)
(*    Enrico Bertolazzi                                                *)
(*    Department of Industrial Engineering                             *)
(*    University of Trento                                             *)
(*    e-mail: enrico.bertolazzi@unitn.it                               *)
(*                                                                     *)
(***********************************************************************)
*/

///
/// file: CollisionDetection_AlignedBox.cc
///

#include "CollisionDetection_AlignedBox.hh"

namespace CollisionDetection
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  AlignedBox<D> &
  AlignedBox<D>::operator=(
    AlignedBox<D> const &AlignedBox_in)
  {
    if (this == &AlignedBox_in)
    {
      return *this;
    }
    else
    {
      this->min() = AlignedBox_in.min();
      this->max() = AlignedBox_in.max();
      return *this;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  AlignedBox<D>::isDegenerated(
    Real tolerance)
    const
  {
    return this->m_min.isApprox(this->m_max, tolerance);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Vector<D>
  AlignedBox<D>::centroid(void)
    const
  {
    return this->center();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  AlignedBox<D>::centerDistance(
    Vector<D> const &point)
    const
  {
    if (this->contains(point))
      return Real(0.0);
    else
      return (this->center() - point).norm();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  AlignedBox<D>::centerDistance(
    AlignedBox<D> const &box)
    const
  {
    if (this->collides(box))
      return Real(0.0);
    else
      return (this->center() - box.center()).norm();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} // namespace CollisionDetection

///
/// eof: CollisionDetection_AlignedBox.cc
///