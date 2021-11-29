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
/// file: CollisionDetection_AlignedBox.hxx
///

#pragma once

#ifndef INCLUDE_COLLISION_DETECTION_ALIGNEDBOX
#define INCLUDE_COLLISION_DETECTION_ALIGNEDBOX

#include "CollisionDetection.hh"

namespace CollisionDetection
{

  /*\
   |      _    _ _                      _ ____
   |     / \  | (_) __ _ _ __   ___  __| | __ )  _____  __
   |    / _ \ | | |/ _` | '_ \ / _ \/ _` |  _ \ / _ \ \/ /
   |   / ___ \| | | (_| | | | |  __/ (_| | |_) | (_) >  <
   |  /_/   \_\_|_|\__, |_| |_|\___|\__,_|____/ \___/_/\_\
   |               |___/
  \*/

  //! Axis-aliged box class container
  /**
   * This class represents an axis-aligned box as a pair of
   * the minimal and maximal corners.
   */
  template <size_t D>
  class AlignedBox : public Eigen::AlignedBox<Real, D>
  {
  public:
    /**
     * Assignment operator=.
     */
    AlignedBox &
    operator=(
      const AlignedBox &AlignedBox_in //!< Input AlignedBox object
    );

    /**
     * \returns \c true if \c *this is degenerated to a point, within the precision
     * determined by \a tolerance.
     */
    bool
    isDegenerated(
      Real tolerance = EPSILON //!< Tolerance
    ) const;

    /**
     * \returns the centroid of the box.
     */
    Vector<D>
    centroid(void) const;

    /**
     * \returns the distance between the point \a point and the box center \c *this,
     * and zero if \a point is inside the box.
     */
    Real
    centerDistance(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns the distance between the boxes \a box and \c *this centers,
     * and zero if the boxes intersect.
     */
    Real
    centerDistance(
      const AlignedBox &box //!< Query AlignedBox object
    ) const;

  }; // class AlignedBox

} // namespace CollisionDetection

#endif // INCLUDE_COLLISION_DETECTION_ALIGNEDBOX

///
/// eof: CollisionDetection_AlignedBox.hh
///