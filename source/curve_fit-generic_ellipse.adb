--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

package body Curve_Fit.Generic_Ellipse is

   function Zero return Number is (One - One);

   function "-" (L : Number) return Number is (Zero - L);

   function "abs" (L : Number) return Number is
      (if L < Zero then -L else L);

   subtype Small_Int is Positive range 2 .. 3;

   function "**" (L : Number; R : Small_Int) return Number;
   function "*" (L : Small_Int; R : Number) return Number;

   function "*" (L : Small_Int; R : Number) return Number is
      R_2 : constant Number := R + R;
   begin
      return (if L = 2 then R_2 else R_2 + R);
   end "*";

   ----------
   -- "**" --
   ----------

   function "**" (L : Number; R : Small_Int) return Number is
      L_2 : constant Number := L * L;
   begin
      return (if R = 2 then L_2 else L_2 * L);
   end "**";

   -----------------
   -- Ellipse_Fit --
   -----------------

   procedure Ellipse_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_List;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50)
   is
      subtype Matrix_T is Matrix (1 .. 5, Points'Range);

      A : Number renames Result (Semi_Major_Axis);
      B : Number renames Result (Semi_Minor_Axis);
      T : Number renames Result (Tilt_Angle);

      function Inside_Ellipse (U, V : Number) return Boolean is
        (U ** 2 / A ** 2 + V ** 2 / B ** 2 < One);

      Frame   : Frame_Parameters renames Result (Frame_Parameter_Index);
      Ellipse : Ellipse_Parameters renames Result (Ellipse_Parameter_Index);
   begin
      Result := Initial;
      RSS := -One;

      for Step in 1 .. Max_Steps loop
         declare
            J_T : Matrix_T;
            M   : Matrix (1 .. 5, 1 .. 5) := [others => [others => Zero]];
            Residuals : Vector (Points'Range);
         begin
            --  Build J_T and Residuals
            for K in Points'Range loop
               declare
                  P  : constant Vector_2D := Points (K);
                  CP : constant Vector_2D := To_Canonical_Frame (P, Frame);

                  UV : constant Vector_2D :=
                    Ellipse_Projection (CP, Ellipse, Epsilon);

                  Pr : constant Vector_2D := From_Canonical_Frame (UV, Frame);

                  U : constant Number := UV (1);
                  V : constant Number := UV (2);

                  Sin_T : constant Number := Sin (T);
                  Cos_T : constant Number := Cos (T);

                  dP_dx : constant Number :=
                    (2 * Cos_T * U / A ** 2 - 2 * Sin_T * V / B ** 2);
                  --  ∂P / ∂x

                  dP_dy : constant Number :=
                    (2 * Sin_T * U / A ** 2 + 2 * Cos_T * V / B ** 2);
                  --  ∂P / ∂y

                  dP_dx0 : constant Number :=
                    (2 * Sin_T * V / B ** 2 - 2 * Cos_T * U / A ** 2);
                  --  ∂P / ∂x₀

                  dP_dy0 : constant Number :=
                    (Zero - 2 * Sin_T * U / A ** 2 - 2 * Cos_T * V / B ** 2);
                  --  ∂P / ∂y₀

                  dP_dt : constant Number :=
                    (2 * U * V / A ** 2 - 2 * U * V / B ** 2);
                  --  ∂P / ∂θ

                  dP_da : constant Number :=
                    (Zero - 2 * U ** 2 / A ** 3);
                  --  ∂P / ∂a

                  dP_db : constant Number :=
                    (Zero - 2 * V ** 2 / B ** 3);
                  --  ∂P / ∂b

                  XY : constant Number := Sqrt (dP_dx ** 2 + dP_dy ** 2);

                  RV : constant Vector_2D := P - Pr;
               begin
                  J_T (1, K) := dP_dx0 / XY;
                  J_T (2, K) := dP_dy0 / XY;
                  J_T (3, K) := dP_dt / XY;
                  J_T (4, K) := dP_da / XY;
                  J_T (5, K) := dP_db / XY;

                  Residuals (K) := Sqrt (RV (1) ** 2 + RV (2) ** 2);

                  if Inside_Ellipse (CP (1), CP (2)) then
                     Residuals (K) := -Residuals (K);
                  end if;
               end;
            end loop;

            for K in Matrix_T'Range (2) loop
               for I in Matrix_T'Range (1) loop
                  for J in Matrix_T'Range (1) loop
                     M (I, J) := @ + J_T (I, K) * J_T (J, K);
                  end loop;
               end loop;
            end loop;

            M := Inverse (M);

            J_T := M * J_T;

            declare
               D : constant Parameters := Parameters (J_T * Residuals);
            begin
               for J in D'Range loop
                  Result (J) := @ - D (J);
               end loop;
            end;

            declare
               Next_RSS : constant Number :=
                 [for R of Residuals => R ** 2]'Reduce ("+", Zero);
            begin
               exit when abs (RSS - Next_RSS) < Epsilon;
               RSS := Next_RSS;
            end;
         end;
      end loop;
   end Ellipse_Fit;

   ------------------------
   -- Ellipse_Projection --
   ------------------------

   function Ellipse_Projection
     (Point      : Vector_2D;
      Parameters : Ellipse_Parameters;
      Epsilon    : Number) return Vector_2D
   is
      CP : Vector_2D := Point;
      U  : Number renames CP (1);
      V  : Number renames CP (2);
      A  : Number renames Parameters (Semi_Major_Axis);
      B  : Number renames Parameters (Semi_Minor_Axis);

      function F (T : Number) return Number;
      function DF (T : Number) return Number;

      function F (T : Number) return Number is
         P1 : constant Number := A * U / (T + A * A);
         P2 : constant Number := B * V / (T + B * B);
      begin
         return P1 ** 2 + P2 ** 2 - One;
      end F;

      function DF (T : Number) return Number is
         P1    : constant Number := 2 * (A * U) ** 2 / (T + A ** 2) ** 3;
         P2    : constant Number := 2 * (B * V) ** 2 / (T + B ** 2) ** 3;
      begin
         return -P1 - P2;
      end DF;

      Turn : array (Point'Range) of Boolean := [False, False];
      T    : Number;
      S    : Number;
   begin
      --  Mirror CP to the first quadrant
      for J in CP'Range loop
         if CP (J) < Zero then
            CP (J) := -CP (J);
            Turn (J) := True;
         end if;
      end loop;

      --  Initialize T: T₀ = max {au-a**2, bv-b**2}
      declare
         T1 : constant Number := A * (U - A);
         T2 : constant Number := B * (V - B);
      begin
         T := (if T1 < T2 then T2 else T1);
      end;

      --  Find T with Newton method
      for J in 1 .. 10 loop
         S := F (T) / DF (T);
         T := T - S;

         exit when abs S < Epsilon;
      end loop;

      --  Calculate the point from T
      CP :=
        [A ** 2 * U / (T + A ** 2),
         B ** 2 * V / (T + B ** 2)];

      --  Mirror CP back
      for J in CP'Range loop
         if Turn (J) then
            CP (J) := -CP (J);
         end if;
      end loop;

      return CP;
   end Ellipse_Projection;

   --------------------------
   -- From_Canonical_Frame --
   --------------------------

   function From_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D
   is
      R : constant Vector_2D :=
        [1 => V (1) * Cos (P (Tilt_Angle)) -
              V (2) * Sin (P (Tilt_Angle)),
         2 => V (1) * Sin (P (Tilt_Angle)) +
              V (2) * Cos (P (Tilt_Angle))];
      C : constant Vector_2D := Center (P);
   begin
      return C + R;
   end From_Canonical_Frame;

   ------------------------
   -- To_Canonical_Frame --
   ------------------------

   function To_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D
   is
      R : constant Vector_2D := V - Center (P);
   begin
      return [1 => R (1) * Cos (P (Tilt_Angle)) +
                   R (2) * Sin (P (Tilt_Angle)),
              2 => R (2) * Cos (P (Tilt_Angle)) -
                   R (1) * Sin (P (Tilt_Angle))];
   end To_Canonical_Frame;

end Curve_Fit.Generic_Ellipse;
