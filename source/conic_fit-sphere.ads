--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Ada.Numerics.Elementary_Functions;
with Ada.Numerics.Real_Arrays;

use Ada.Numerics.Elementary_Functions;
use Ada.Numerics.Real_Arrays;

with Conic_Fit.Float_Vectors;
with Conic_Fit.Generic_Sphere;

package Conic_Fit.Sphere is new Conic_Fit.Generic_Sphere
  (Number  => Float,
   Zero    => 0.0,
   One     => 1.0,
   Matrix  => Real_Matrix,
   Vector  => Real_Vector,
   Vectors => Conic_Fit.Float_Vectors);

pragma Pure (Conic_Fit.Sphere);
