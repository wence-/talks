Point(1) = {1.0, 0, 0, 0.1};
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{5};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{6};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{7};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{8};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{9};
}
Point(10) = {100.0, 0, 0, 2.0};
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{10};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{11};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{12};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{13};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{14};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{15};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{16};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{17};
}
Line Loop(17) = {16, 9, 10, 11, 12, 13, 14, 15};
Line Loop(18) = {8, 1, 2, 3, 4, 5, 6, 7};
Plane Surface(19) = {17, 18};
Physical Surface(20) = {19};
Physical Line(21) = {2, 1, 8, 7, 6, 5, 4, 3};
Physical Line(22) = {10, 9, 16, 15, 14, 13, 12, 11};
