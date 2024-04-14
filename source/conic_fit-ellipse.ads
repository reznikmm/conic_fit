--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Ada.Numerics.Elementary_Functions;
with Ada.Numerics.Real_Arrays;

use Ada.Numerics.Elementary_Functions;
use Ada.Numerics.Real_Arrays;

with Conic_Fit.Generic_Ellipse;

package Conic_Fit.Ellipse is new Conic_Fit.Generic_Ellipse
  (Number => Float,
   One    => 1.0,
   Matrix => Real_Matrix,
   Vector => Real_Vector);

pragma Pure (Conic_Fit.Ellipse);
