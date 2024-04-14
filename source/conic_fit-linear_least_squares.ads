--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

--  This is a generic version of the linear least squares solver.
--
--  For given A*X = B (where X, B - vectors and A is a matrix such as
--  A'Length (1) > A'Length (2)), it finds X by minimizing ||A*X - B||**2
--
--  We use functions that provide A and B values to reduce memory consumption.
--
--  Vector/matrix element is a private type (not a float) to be able to use
--  fixed point types.
--

generic
   type Element is private;
   type Element_Vector is array (Integer range <>) of Element;

   type Element_Matrix is
     array (Integer range <>, Integer range <>) of Element;

   with function Inverse (M : Element_Matrix) return Element_Matrix is <>;
   --  Do matrix inversion. M'Range (1) = M'Range (2) = X'Range

   with function A (J, K : Positive) return Element;
   --  J in 1 .. Rows, K in X'Range

   with function B (J : Positive) return Element;
   --  J in 1 .. Rows

   with function "*" (L, R : Element) return Element is <>;
   with function "+" (L, R : Element) return Element is <>;

   Zero : Element;
procedure Conic_Fit.Linear_Least_Squares
  (Rows : Positive;
   X    : out Element_Vector);
--  Solve A*X = B, when A and B are given as functions

pragma Pure (Conic_Fit.Linear_Least_Squares);
