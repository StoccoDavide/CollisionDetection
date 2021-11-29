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
/// file: CollisionDetection_Simplex.cc
///

#include "CollisionDetection_Simplex.hh"

namespace CollisionDetection
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D>::~Simplex(void)
  {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D>::Simplex(void)
  {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D>::Simplex(
    Matrix<D, D + 1> const &points) : m_points(points)
  {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D> &
  Simplex<D>::operator=(
    Simplex<D> const &other)
  {
    if (this == &other)
    {
      return *this;
    }
    else
    {
      this->m_points = other.m_points;
      return *this;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  size_t
  Simplex<D>::dim(void)
    const
  {
    return size_t(this->m_points.rows());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  size_t
  Simplex<D>::size(void)
    const
  {
    return size_t(this->m_points.cols());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::isEmpty(void)
    const
  {
    return this->m_points.array().any();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  void
  Simplex<D>::setEmpty(void)
  {
    this->m_points.setConstant(QUIET_NAN);
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Vector<D> const &
  Simplex<D>::point(
    size_t i)
    const
  {
    return this->m_points.col(i);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Vector<D> &
  Simplex<D>::point(
    size_t i)
  {
    return this->m_points.col(i);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Vector<D>
  Simplex<D>::centroid(void)
    const
  {
    return this->m_points.colwise().sum() / Real(D);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Vector<D>
  Simplex<D>::center(void)
    const
  {
    return this->center();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Vector<D + 1>
  Simplex<D>::barycentric(
    const Vector<D> &point)
    const
  {
    Vector<D>     lambda((this->m_points.block(0, 0, D, D) - this->m_points.cols(Eigen::last)).inverse() * (point - this->m_points.cols(Eigen::last)));
    Vector<D + 1> result;
    result << lambda, Real(1.0) - lambda.rowwise().sum();
    return result;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::contains(
    const Vector<D> &point)
    const
  {
    Vector<D + 1> result(this->barycentric(point));
    return (result.array() >= Real(0.0)).any() && result.rowwise().sum() <= Real(D);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::contains(
    const Simplex<D> &simplex)
    const
  {
    for (size_t i = 0; i < D; ++i)
    {
      if (!this->contains(simplex.point(i)))
        return false;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::intersects(
    const Simplex<D> &simplex)
    const
  {
    for (size_t i = 0; i < D; ++i)
    {
      if (this->contains(simplex.point(i)))
        return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D> &
  Simplex<D>::translate(
    const Vector<D> &translation)
  {
    this->m_points.colwise() += translation;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Simplex<D>
  Simplex<D>::translated(
    const Vector<D> &translation)
    const
  {
    Simplex<D> result(this->m_points);
    result.translate(translation);
    return result;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  template <int Mode, int Options>
  void
  Simplex<D>::transform(
    const Eigen::Transform<Real, D, Mode, Options> &transform)
  {
    this->translate(transform.translation());
    this->m_points.colwise() * transform.linear();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  template <int Mode, int Options>
  Simplex<D>
  Simplex<D>::transformed(
    const Eigen::Transform<Real, D, Mode, Options> &transform)
    const
  {
    Simplex<D> result(this->m_points);
    result.transform(transform);
    return result;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::isApprox(
    Simplex<D> const &other,
    Real              tolerance)
    const
  {
    return this->m_points.isApprox(other.m_points);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  bool
  Simplex<D>::isDegenerated(
    Real tolerance)
    const
  {
    return this->m_point.col(0).isApprox(this->center(), tolerance);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::exteriorDistance(
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
  Simplex<D>::exteriorDistance(
    Simplex<D> const &simplex)
    const
  {
    if (this->collides(simplex))
      return Real(0.0);
    else
      return (this->center() - simplex.center()).norm();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::squaredExteriorDistance(
    const Vector<D> &point)
    const
  {
    if (this->contains(point))
      return Real(0.0);
    else
      return Real(0.0); // TODO
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::squaredExteriorDistance(
    const Simplex &simplex)
    const
  {
    if (this->collides(simplex))
      return Real(0.0);
    else
      return Real(0.0); // TODO
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::exteriorDistance(
    const Vector<D> &point)
    const
  {
    if (this->contains(point))
      return Real(0.0);
    else
      return Real(0.0); // TODO
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::exteriorDistance(
    const Simplex &simplex)
    const
  {
    if (this->collides(simplex))
      return Real(0.0);
    else
      return Real(0.0); // TODO
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t D>
  Real
  Simplex<D>::centerDistance(
    const Vector<D> &point)
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
  Simplex<D>::centerDistance(
    const Simplex &simplex)
    const
  {
    if (this->collides(simplex))
      return Real(0.0);
    else
      return (this->center() - simplex.center()).norm();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} // namespace CollisionDetection

///
/// eof: CollisionDetection_Simplex.cc
///