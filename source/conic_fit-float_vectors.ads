--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Ada.Numerics.Real_Arrays;

with Conic_Fit.Generic_Vectors;

package Conic_Fit.Float_Vectors is new Conic_Fit.Generic_Vectors
  (Float, Ada.Numerics.Real_Arrays.Real_Vector);

pragma Pure (Conic_Fit.Float_Vectors);
