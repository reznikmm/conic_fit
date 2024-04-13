--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Ada.Float_Text_IO;
with Ada.Integer_Text_IO;
with Ada.Text_IO;
with Ada.Containers.Vectors;

package body PLY is
   use type Point;

   package Point_Vectors is new Ada.Containers.Vectors (Positive, Point);

   type Indexes is array (1 .. 3) of Positive;

   package Face_Vectors is new Ada.Containers.Vectors (Positive, Indexes);

   function To_Index (P : Point) return Positive;

   Vertex : Point_Vectors.Vector;
   Faces  : Face_Vectors.Vector;

   function To_Index (P : Point) return Positive is
      X : Natural := Vertex.Find_Index (P);
   begin
      if X = Point_Vectors.No_Index then
         Vertex.Append (P);
         X := Vertex.Last_Index;
      end if;

      return X;
   end To_Index;

   procedure Append_Face (F : Face) is
      V : constant Indexes := [for P of F => To_Index (P)];
   begin
      Faces.Append (V);
   end Append_Face;

   procedure Write
     (Scale : Float := 1.0;
      Name  : String := "/tmp/aaa.ply")
   is
      use Ada.Text_IO;
      Output : File_Type;
   begin
      Create (Output, Name => Name);
      Put_Line (Output, "ply");
      Put_Line (Output, "format ascii 1.0");
      Put (Output, "element vertex");
      Ada.Integer_Text_IO.Put (Output, Vertex.Last_Index);
      New_Line (Output);
      Put_Line (Output, "property float x");
      Put_Line (Output, "property float y");
      Put_Line (Output, "property float z");
      Put (Output, "element face");
      Ada.Integer_Text_IO.Put (Output, Faces.Last_Index);
      New_Line (Output);
      Put_Line (Output, "property list uint uint vertex_indices");
      Put_Line (Output, "end_header");

      for V of Vertex loop
         Ada.Float_Text_IO.Put (Output, V (1) * Scale, Aft => 8, Exp => 0);
         Put (Output,  ' ');
         Ada.Float_Text_IO.Put (Output, V (2) * Scale, Aft => 8, Exp => 0);
         Put (Output,  ' ');
         Ada.Float_Text_IO.Put (Output, V (3) * Scale, Aft => 8, Exp => 0);
         New_Line (Output);
      end loop;

      for F of Faces loop
         Ada.Integer_Text_IO.Put (Output, 3);
         Ada.Integer_Text_IO.Put (Output, F (1) - 1);
         Ada.Integer_Text_IO.Put (Output, F (2) - 1);
         Ada.Integer_Text_IO.Put (Output, F (3) - 1);
         New_Line (Output);
      end loop;

      Close (Output);
      Faces.Clear;
   end Write;
end PLY;
