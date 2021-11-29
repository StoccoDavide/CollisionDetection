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
/// file: CollisionDetection_Types.hh
///

#pragma once

#ifndef INCLUDE_COLLISION_DETECTION_TYPES
#define INCLUDE_COLLISION_DETECTION_TYPES

#include "CollisionDetection.hh"

namespace CollisionDetection
{
  /*\
   |   _____
   |  |_   _|   _ _ __   ___  ___
   |    | || | | | '_ \ / _ \/ __|
   |    | || |_| | |_) |  __/\__ \
   |    |_| \__, | .__/ \___||___/
   |        |___/|_|
  \*/

  typedef double       Real;       //!< Real number type
  typedef int          integer;    //!< Integer number type
  typedef std::ostream out_stream; //!< Output stream type

  //! Nx1 matrix of real number alias
  template <size_t M>
  using Vector = Eigen::Matrix<Real, M, 1>;

  //! MxN matrix of real number alias
  template <size_t M, size_t N>
  using Matrix = Eigen::Matrix<Real, M, N>;

  typedef Eigen::Matrix<Real, 2, 1>                           vec2; //!< 2x1 vector type (column vector)
  typedef Eigen::Matrix<Real, 2, 2>                           mat2; //!< 2x2 matrix type
  typedef Eigen::Matrix<Real, 3, 1>                           vec3; //!< 3x1 vector type (column vector)
  typedef Eigen::Matrix<Real, 3, 3>                           mat3; //!< 3x3 matrix type
  typedef Eigen::Matrix<Real, 4, 1>                           vec4; //!< 4x1 vector type (column vector)
  typedef Eigen::Matrix<Real, 4, 4>                           mat4; //!< 4x4 matrix type
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>              vecN; //!< Nx1 vector of Real number type (column vector)
  typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> matN; //!< NxN matrix of Real number type

  static Real const EPSILON_MACHINE = std::numeric_limits<Real>::epsilon();      //!< Machine epsilon epsilon static constant value
  static Real const EPSILON_HIGH    = 1.0E-16;                                   //!< High precision epsilon static constant value
  static Real const EPSILON_MEDIUM  = 1.0E-10;                                   //!< Medium precision epsilon static constant value
  static Real const EPSILON_LOW     = 1.0E-07;                                   //!< Low precision epsilon static constant value
  static Real const EPSILON         = EPSILON_MEDIUM;                            //!< Standard precision epsilon static constant value
  static Real const INFTY           = std::numeric_limits<Real>::infinity();     //!< Infinity static constant value
  static Real const QUIET_NAN       = std::numeric_limits<Real>::quiet_NaN();    //!< Not-a-Number static constant value
  static Real const PI              = Real(3.141592653589793238462643383279500); //!< Pi static constant value
  static Real const PIDIV180        = Real(0.017453292519943295769236907684886); //!< Pi/180 static constant value

} // namespace CollisionDetection

#endif // INCLUDE_COLLISION_DETECTION_TYPES

///
/// file: CollisionDetection_Types.hh
///