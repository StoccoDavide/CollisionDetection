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
/// file: CollisionDetection_Simplex.hh
///

#pragma once

#ifndef INCLUDE_COLLISION_DETECTION_SIMPLEX
#define INCLUDE_COLLISION_DETECTION_SIMPLEX

#include "CollisionDetection.hh"

namespace CollisionDetection
{

  /*\
   |   ____  _                 _
   |  / ___|(_)_ __ ___  _ __ | | _____  __
   |  \___ \| | '_ ` _ \| '_ \| |/ _ \ \/ /
   |   ___) | | | | | | | |_) | |  __/>  <
   |  |____/|_|_| |_| |_| .__/|_|\___/_/\_\
   |                    |_|
  \*/

  //! N-dimensional simplex class container
  /**
   * This class represents a simplex defined by a number of points N
   * and by the dimensions D of the space in which the simplex is
   * defined.
   */
  template <size_t D>
  class Simplex
  {
  private:
    Matrix<D, D + 1> m_points; //!< Simplex points

    /**
     * Simplex class move constructor.
     */
    Simplex(Simplex &&) = delete;

  public:
    /**
     * Simplex class copy constructor.
     */
    Simplex(const Simplex &) = default;

    /**
     * Simplex class destructor.
     */
    ~Simplex(void);

    /**
     * Simplex class constructor.
     */
    Simplex(void);

    /**
     * Simplex class constructor.
     */
    Simplex(
      Matrix<D, D + 1> const &points //!< Input points
    );

    /**
     * Assignment operator=.
     */
    Simplex &
    operator=(
      Simplex const &other //!< Input Simplex object
    );

    /**
     * \returns the dimension in which the simplex holds.
     */
    size_t
    dim(void) const;

    /**
     * \returns the simplex points number.
     */
    size_t
    size(void) const;

    /**
     * \returns true if the simplex is empty.
     */
    bool
    isEmpty(void) const;

    /**
     * Makes \c *this an empty simplex.
     */
    void
    setEmpty(void);

    /**
     * \returns simplex i-th point const reference.
     */
    Vector<D> const &
    point(
      size_t i //!< Vector index
    ) const;

    /**
     * \returns simplex i-th point reference.
     */
    Vector<D> &
    point(
      size_t i //!< Vector index
    );

    /**
     * \returns the centroid of the simplex.
     */
    Vector<D>
    centroid(void) const;

    /**
     * \returns the center (centroid) of the simplex.
     */
    Vector<D>
    center(void) const;

    /**
     * \returns the barycentric coordinates of the point \a point.
     */
    Vector<D + 1>
    barycentric(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns true if the point \a point is inside the simplex \c *this.
     */
    bool
    contains(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns true if the simplex \a simplex is entirely inside the simplex \c *this.
     */
    bool
    contains(
      const Simplex &simplex //! Input simplex object
    ) const;

    /**
     * \returns true if the simplex \a simplex is intersecting the simplex \c *this.
     */
    bool
    intersects(
      const Simplex &simplex //! Input simplex object
    ) const;

    /**
     * Translate \c *this by the vector \a translation and returns a reference to \c *this.
     */
    Simplex &
    translate(
      const Vector<D> &translation //!< Translation vector
    );

    /**
     * \returns a copy of \c *this translated by the vector \a translation.
     */
    Simplex
    translated(
      const Vector<D> &translation //!< Translation vector
    ) const;

    /**
     * Transformed \c *this by transformation \a transform.
     */
    template <int Mode, int Options>
    void
    transform(
      const Eigen::Transform<Real, D, Mode, Options> &transform //!< Input Transformation object
    );

    /**
     * \returns a copy of \c *this transformed by a transformation \a transform.
     */
    template <int Mode, int Options>
    Simplex
    transformed(
      const Eigen::Transform<Real, D, Mode, Options> &transform //!< Input Transformation object
    ) const;

    /**
     * \returns \c true if \c *this is approximately equal to \a other, within the precision
     * determined by \a tolerance.
     */
    bool
    isApprox(
      const Simplex &other,              //!< Input Simplex object
      Real           tolerance = EPSILON //!< Tolerance
    ) const;

    /**
     * \returns \c true if \c *this is degenerated to a point, within the precision
     * determined by \a tolerance.
     */
    bool
    isDegenerated(
      Real tolerance = EPSILON //!< Tolerance
    ) const;

    /**
     * \returns the squared distance between the point \a point and the simplex \c *this,
     * and zero if \a point is inside the simplex.
     */
    Real
    squaredExteriorDistance(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns the squared distance between the simplexes \a simplex and \c *this,
     * and zero if the simplexes intersect.
     */
    Real
    squaredExteriorDistance(
      const Simplex &simplex //!< Query simplex object
    ) const;

    /**
     * \returns the distance between the simplexes \a simplex and \c *this,
     * and zero if the simplexes intersect.
     */
    Real
    exteriorDistance(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns the distance between the simplexes \a simplex and \c *this,
     * and zero if the simplexes intersect.
     */
    Real
    exteriorDistance(
      const Simplex &simplex //!< Query simplex object
    ) const;

    /**
     * \returns the distance between the point \a point and the simplex center \c *this,
     * and zero if \a point is inside the simplex.
     */
    Real
    centerDistance(
      const Vector<D> &point //!< Query point
    ) const;

    /**
     * \returns the distance between the simplexes \a simplex and \c *this centers,
     * and zero if the simplexes intersect.
     */
    Real
    centerDistance(
      const Simplex &simplex //!< Query simplex object
    ) const;

  }; // class Simplex

} // namespace CollisionDetection

#endif // INCLUDE_COLLISION_DETECTION_SIMPLEX

///
/// eof: CollisionDetection_Simplex.hh
///