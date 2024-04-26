--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Ada.Unchecked_Deallocation;

package body Conic_Fit.Generic_Ellipsoid is

   function Zero return Number is (One - One);
   function Two return Number is (One + One);

   function "-" (L : Number) return Number is (Zero - L);

   function "abs" (L : Number) return Number is
      (if L < Zero then -L else L);

   subtype Small_Int is Positive range 2 .. 3;

   function "**" (L : Number; R : Small_Int) return Number;
   function "*" (L : Small_Int; R : Number) return Number;

   function "+" (L, R : Vector_3D) return Vector_3D is
     [L (1) + R (1), L (2) + R (2), L (3) + R (3)];

   function "-" (L, R : Vector_3D) return Vector_3D is
     [L (1) - R (1), L (2) - R (2), L (3) - R (3)];

   function Center (P : Frame_Parameters) return Vector_3D is
     [P (Center_X), P (Center_Y), P (Center_Z)];

   ---------
   -- "*" --
   ---------

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
   -- Ellipsoid_Fit --
   -----------------

   procedure Ellipsoid_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_Array;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50)
   is
      subtype Matrix_T is Matrix (1 .. 9, Points'Range);
      --  A type of the Jacobian matrix (transposed).

      subtype Residual_Vector is Vector (Points'Range);

      type Heap_Item is record
         J_T       : Matrix_T;
         Residuals : Residual_Vector;
      end record;

      type Heap_Item_Access is access Heap_Item;

      procedure Free is new Ada.Unchecked_Deallocation
        (Heap_Item, Heap_Item_Access);

      subtype Mx9x9 is Matrix (Matrix_T'Range (1), Matrix_T'Range (1));
      subtype Vx9 is Vector (Matrix_T'Range (1));

      procedure MxJ (M : Mx9x9; J : in out Matrix_T);
      --  Calculate J := M * J;

      function "*" (J_T : Matrix_T; V : Vector) return Vx9;
      function Find_RSS (Copy : Parameters) return Number;

      ---------
      -- "*" --
      ---------

      function "*" (J_T : Matrix_T; V : Vector) return Vx9 is
         function Cell (R : Positive) return Number;

         function Cell (R : Positive) return Number is
           ([for C in V'Range => J_T (R, C) * V (C)]'Reduce ("+", Zero));
      begin
         return [for R in Matrix_T'Range (1) => Cell (R)];
      end "*";

      --------------
      -- Find_RSS --
      --------------

      function Find_RSS (Copy : Parameters) return Number is
         Frame     : Frame_Parameters renames Copy (Frame_Parameters'Range);
         Ellipsoid : Ellipsoid_Parameters renames
           Copy (Ellipsoid_Parameter_Index);

         Total : Number := Zero;
      begin
         --  Build Residuals
         for K in Points'Range loop
            declare
               P  : constant Vector_3D := Points (K);
               CP : constant Vector_3D := To_Canonical_Frame (P, Frame);

               UVW : constant Vector_3D :=
                 Ellipsoid_Projection (CP, Ellipsoid, Epsilon);

               Pr : constant Vector_3D := From_Canonical_Frame (UVW, Frame);
               RV : constant Vector_3D := P - Pr;
               --  Residual vector

               Residual : Number;
            begin
               Residual := Sqrt (RV (1) ** 2 + RV (2) ** 2 + RV (3) ** 2);
               Total := @ + Residual ** 2;
            end;
         end loop;

         return Total;
      end Find_RSS;

      procedure MxJ (M : Mx9x9; J : in out Matrix_T) is
      begin
         for C in J'Range (2) loop  --  column J [1 .. Points'Length]
            declare
               V : constant Vector (J'Range (1)) :=  --  Column J(C) copy
                 [for K in J'Range (1) => J (K, C)];
            begin
               for R in J'Range (1) loop  --  row J [1 .. 9]
                  J (R, C) :=
                    [for K in J'Range (1) => M (R, K) * V (K)]
                      'Reduce ("+", Zero);
               end loop;
            end;
         end loop;
      end MxJ;

      Heap      : Heap_Item_Access := new Heap_Item;
      Frame     : Frame_Parameters renames Result (Frame_Parameters'Range);
      Ellipsoid : Ellipsoid_Parameters renames
        Result (Ellipsoid_Parameter_Index);

   begin
      Result := Initial;
      RSS := -One;

      Top_Loop :
      for Step in 1 .. Max_Steps loop
         declare
            A     : constant Number := Result (Semi_Major_Axis);
            B     : constant Number := Result (Semi_Middle_Axis);
            C     : constant Number := Result (Semi_Minor_Axis);
            Sin_A : constant Number := Sin (Result (Roll));
            Cos_A : constant Number := Cos (Result (Roll));
            Sin_B : constant Number := Sin (Result (Pitch));
            Cos_B : constant Number := Cos (Result (Pitch));
            Sin_C : constant Number := Sin (Result (Yaw));
            Cos_C : constant Number := Cos (Result (Yaw));

            function Inside_Ellipsoid (U, V, W : Number) return Boolean is
              ((U / A) ** 2 + (V / B) ** 2 + (W / C) ** 2 < One);
            --  u²/a² + v²/b² + w²/c²  < 1

            J_T : Matrix_T renames Heap.J_T;
            M   : Matrix (Matrix_T'Range (1), Matrix_T'Range (1)) :=
              [others => [others => Zero]];

            Residuals : Residual_Vector renames Heap.Residuals;
            Change    : Vector (Matrix_T'Range (1));
            Scale     : Number := One;
         begin
            --  Build J_T and Residuals
            for K in Points'Range loop
               declare
                  P  : constant Vector_3D := Points (K);
                  CP : constant Vector_3D := To_Canonical_Frame (P, Frame);

                  UVW : constant Vector_3D :=
                    Ellipsoid_Projection (CP, Ellipsoid, Epsilon);

                  Pr : constant Vector_3D := From_Canonical_Frame (UVW, Frame);

                  CX : constant Number := Pr (1) - Result (Center_X);
                  CY : constant Number := Pr (2) - Result (Center_Y);
                  CZ : constant Number := Pr (3) - Result (Center_Z);

                  U : constant Number := UVW (1);
                  V : constant Number := UVW (2);
                  W : constant Number := UVW (3);

                  dP_dx : constant Number :=
                    2 * U * Cos_B * Cos_C / A ** 2 +
                    2 * V * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C) / B ** 2 +
                    2 * W * (Sin_A * Sin_C - Cos_A * Sin_B * Cos_C) / C ** 2;
                  --  ∂P / ∂x

                  dP_dy : constant Number :=
                    2 * U * (-Cos_B * Sin_C) / A ** 2 +
                    2 * V * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C) / B ** 2 +
                    2 * W * (Sin_A * Cos_C + Cos_A * Sin_B * Sin_C) / C ** 2;
                  --  ∂P / ∂y

                  dP_dz : constant Number :=
                    2 * U * Sin_B / A ** 2 +
                    2 * V * (-Sin_A * Cos_B) / B ** 2 +
                    2 * W * Cos_A * Cos_B / C ** 2;
                  --  ∂P / ∂z

                  XYZ   : constant Number :=
                    Sqrt (dP_dx ** 2 + dP_dy ** 2 + dP_dz ** 2);
                  --  Sqrt (dP_dx² + dP_dy² + dP_dz²)

                  dP_dx0 : constant Number := -dP_dx;
                  --  ∂P / ∂x₀

                  dP_dy0 : constant Number := -dP_dy;
                  --  ∂P / ∂y₀

                  dP_dz0 : constant Number := -dP_dz;
                  --  ∂P / ∂z₀

                  dP_droll : constant Number :=
                    2 * V / B ** 2 *
                      (CX * (Cos_A * Sin_B * Cos_C - Sin_A * Sin_C)
                       - CY * (Sin_A *  Cos_C + Cos_A * Sin_B * Sin_C)
                       - CZ * Cos_A * Cos_B) +
                    2 * W / C ** 2 *
                      (CX * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C)
                       + CY * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C)
                       - CZ * Sin_A * Cos_B);
                  --  ∂P / ∂α

                  dP_dpitch : constant Number :=
                    2 * U / A ** 2 *
                      (-CX * Sin_B * Cos_C
                       + CY * Sin_B * Sin_C
                       + CZ * Cos_B) +
                    2 * V / B ** 2 *
                      (CX * (Sin_A * Cos_B * Cos_C)
                       - CY * Sin_A * Cos_B * Sin_C
                       + CZ * Sin_A * Sin_B) +
                    2 * W / C ** 2 *
                      (-CX * Cos_A * Cos_B * Cos_C
                       + CY * Cos_A * Cos_B * Sin_C
                       - CZ * Cos_A * Sin_B);
                  --  ∂P / ∂β

                  dP_dyaw : constant Number :=
                    2 * U / A ** 2 *
                      (-CX * Cos_B * Sin_C
                       - CY * Cos_B * Cos_C) +
                    2 * V / B ** 2 *
                      (CX * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C)
                       - CY * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C)) +
                    2 * W / C ** 2 *
                      (CX * (Sin_A * Cos_C + Cos_A * Sin_B * Sin_C)
                       - CY * (Sin_A * Sin_C - Cos_A * Sin_B * Cos_C));
                  --  ∂P / ∂γ

                  dP_da : constant Number := -(2 * U ** 2 / A ** 3);
                  --  ∂P / ∂a

                  dP_db : constant Number := -(2 * V ** 2 / B ** 3);
                  --  ∂P / ∂b

                  dP_dc : constant Number := -(2 * W ** 2 / C ** 3);
                  --  ∂P / ∂c

                  RV    : constant Vector_3D := P - Pr;
                  --  Residual vector
               begin
                  J_T (1, K) := dP_dx0 / XYZ;
                  J_T (2, K) := dP_dy0 / XYZ;
                  J_T (3, K) := dP_dz0 / XYZ;
                  J_T (4, K) := dP_droll / XYZ;
                  J_T (5, K) := dP_dpitch / XYZ;
                  J_T (6, K) := dP_dyaw / XYZ;
                  J_T (7, K) := dP_da / XYZ;
                  J_T (8, K) := dP_db / XYZ;
                  J_T (9, K) := dP_dc / XYZ;

                  Residuals (K) :=
                    Sqrt (RV (1) ** 2 + RV (2) ** 2 + RV (3) ** 2);

                  if Inside_Ellipsoid (CP (1), CP (2), CP (3)) then
                     Residuals (K) := -Residuals (K);
                  end if;
               end;

               --  Fill upper half of M = J x Jᵀ
               for I in Matrix_T'Range (1) loop
                  for J in I .. Matrix_T'Last (1) loop
                     M (I, J) := @ + J_T (I, K) * J_T (J, K);
                  end loop;
               end loop;
            end loop;

            --  Copy upper half of M to lower half
            for I in Matrix_T'Range (1) loop
               for J in I + 1 .. Matrix_T'Last (1) loop
                  M (J, I) := M (I, J);
               end loop;
            end loop;

            M := Inverse (M);

            MxJ (M, J_T);  --  Calculate J_T := M * J_T;

            if Step = 1 then
               RSS := [for R of Residuals => R ** 2]'Reduce ("+", Zero);
            end if;

            Change := J_T * Residuals;

            for J in reverse 1 .. 15 loop
               declare
                  function "+" (J : Positive)
                    return Ellipsoid_Geometric_Parameter_Index
                      is (Ellipsoid_Geometric_Parameter_Index'Val (J - 1));

                  Copy     : Parameters := Result;
                  Next_RSS : Number;
               begin
                  for J in Change'Range loop
                     Copy (+J) := @ - Scale * Change (J);
                  end loop;

                  --  Check ellipsoid parameters invariant
                  if Zero < Copy (Semi_Major_Axis) and then
                    Copy (Semi_Middle_Axis) < Copy (Semi_Major_Axis) and then
                    Copy (Semi_Minor_Axis) < Copy (Semi_Middle_Axis) and then
                    abs (Result (Roll) - Copy (Roll)) < One and then
                    abs (Result (Pitch) - Copy (Pitch)) < One and then
                    abs (Result (Yaw) - Copy (Yaw)) < One
                  then
                     Next_RSS := Find_RSS (Copy);

                     if Next_RSS < RSS then  --  Found better parameters
                        Result := Copy;

                        if abs (RSS - Next_RSS) < Epsilon then
                           RSS := Next_RSS;
                           exit Top_Loop;
                        end if;

                        RSS := Next_RSS;
                        exit;
                     end if;
                  end if;

                  Scale := Scale / Two;  --  Reduce the parameter change
               end;
            end loop;
         end;
      end loop Top_Loop;

      Free (Heap);
   end Ellipsoid_Fit;

   --------------------------
   -- Ellipsoid_Projection --
   --------------------------

   function Ellipsoid_Projection
     (Point      : Vector_3D;
      Parameters : Ellipsoid_Parameters;
      Epsilon    : Number) return Vector_3D
   is
      CP : Vector_3D := Point;
      U  : Number renames CP (1);
      V  : Number renames CP (2);
      W  : Number renames CP (3);
      A  : Number renames Parameters (Semi_Major_Axis);
      B  : Number renames Parameters (Semi_Middle_Axis);
      C  : Number renames Parameters (Semi_Minor_Axis);

      function F (T : Number) return Number;
      function DF (T : Number) return Number;

      function F (T : Number) return Number is
         P1 : constant Number := A * U / (T + A * A);
         P2 : constant Number := B * V / (T + B * B);
         P3 : constant Number := C * W / (T + C * C);
      begin
         return P1 ** 2 + P2 ** 2 + P3 ** 2 - One;
      end F;

      function DF (T : Number) return Number is
         P1 : constant Number := 2 * (A * U) ** 2 / (T + A ** 2) ** 3;
         P2 : constant Number := 2 * (B * V) ** 2 / (T + B ** 2) ** 3;
         P3 : constant Number := 2 * (C * W) ** 2 / (T + C ** 2) ** 3;
      begin
         return -P1 - P2 - P3;
      end DF;

      Turn : array (Point'Range) of Boolean := [others => False];
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

      --  Initialize T: T₀ = max {au-a**2, bv-b**2, cw-c**2}
      declare
         T1 : constant Number := A * (U - A);
         T2 : constant Number := B * (V - B);
         T3 : constant Number := C * (W - C);
      begin
         T := (if T1 < T2 then T2 else T1);
         T := (if T < T3 then T3 else T);
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
         B ** 2 * V / (T + B ** 2),
         C ** 2 * W / (T + C ** 2)];

      --  Mirror CP back
      for J in CP'Range loop
         if Turn (J) then
            CP (J) := -CP (J);
         end if;
      end loop;

      return CP;
   end Ellipsoid_Projection;

   --------------------------
   -- From_Canonical_Frame --
   --------------------------

   function From_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D
   is
      R : constant Vector_3D :=
           [1 => V (1) * Cos (P (Pitch)) * Cos (P (Yaw))
               + V (2) * (Cos (P (Roll)) * Sin (P (Yaw))
                          + Sin (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
               + V (3) * (Sin (P (Roll)) * Sin (P (Yaw))
                          - Cos (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw))),
            2 => -V (1) * Cos (P (Pitch)) * Sin (P (Yaw))
               + V (2) * (Cos (P (Roll)) * Cos (P (Yaw))
                          - Sin (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
               + V (3) * (Sin (P (Roll)) * Cos (P (Yaw))
                          + Cos (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw))),
            3 => V (1) * Sin (P (Pitch))
               - V (2) * Sin (P (Roll)) * Cos (P (Pitch))
               + V (3) * Cos (P (Roll)) * Cos (P (Pitch))];
      C : constant Vector_3D := Center (P);
   begin
      return C + R;
   end From_Canonical_Frame;

   ------------------------
   -- To_Canonical_Frame --
   ------------------------

   function To_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D
   is
      R : constant Vector_3D := V - Center (P);
   begin
      return [1 => R (1) * Cos (P (Pitch)) * Cos (P (Yaw))
                 - R (2) * Cos (P (Pitch)) * Sin (P (Yaw))
                 + R (3) * Sin (P (Pitch)),
              2 => R (1) * (Cos (P (Roll)) * Sin (P (Yaw))
                          + Sin (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
                 + R (2) * (Cos (P (Roll)) * Cos (P (Yaw))
                          - Sin (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
                 - R (3) * Sin (P (Roll)) * Cos (P (Pitch)),
              3 => R (1) * (Sin (P (Roll)) * Sin (P (Yaw))
                          - Cos (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
                 + R (2) * (Sin (P (Roll)) * Cos (P (Yaw))
                          + Cos (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
                 + R (3) * Cos (P (Roll)) * Cos (P (Pitch))];
   end To_Canonical_Frame;

end Conic_Fit.Generic_Ellipsoid;
