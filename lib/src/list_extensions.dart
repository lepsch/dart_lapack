import 'dart:math';

import 'package:lapack/src/matrix.dart';

extension IntFormatExtension on List<int> {
  String get i1 => toString().padLeft(1).w(1);
  String get i2 => toString().padLeft(2).w(2);
  String get i3 => toString().padLeft(3).w(3);
  String get i4 => toString().padLeft(4).w(4);
  String get i5 => toString().padLeft(5).w(5);
  String get i6 => toString().padLeft(6).w(6);
  String get i7 => toString().padLeft(7).w(7);
  String get i8 => toString().padLeft(8).w(8);
  String get i12 => toString().padLeft(12).w(12);
  String get i15 => toString().padLeft(15).w(15);
  String get i36 => toString().padLeft(36).w(36);
}
