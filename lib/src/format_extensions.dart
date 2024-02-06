import 'package:lapack/src/matrix.dart';

extension IntFormatExtension on int {
  String get i2 => toString().padLeft(2);
  String get i3 => toString().padLeft(3);
  String get i4 => toString().padLeft(4);
  String get i5 => toString().padLeft(5);
  String get i6 => toString().padLeft(6);
  String get i7 => toString().padLeft(7);
  String get i8 => toString().padLeft(8);
  String get i12 => toString().padLeft(12);
  String get i36 => toString().padLeft(36);
}

extension DoubleFormatExtension on double {
  String get f4_2 => toStringAsFixed(2).padLeft(4);
  String get f8_2 => toStringAsFixed(2).padLeft(8);

  String get d12_3 => toStringAsFixed(3).padLeft(12);
  String get d12_4 => toStringAsFixed(4).padLeft(12);
  String get d16_6 => toStringAsFixed(6).padLeft(16);
  String get d36_8 => toStringAsFixed(8).padLeft(36);

  String get e15_8 => toStringAsFixed(8).padLeft(15);

  String get g10_3 => toStringAsFixed(3).padLeft(10);
  String get g11_4 => toStringAsFixed(4).padLeft(11);
}

extension StringFormatExtension on String {
  String get a3 => substring(0, 3).padLeft(3);
  String get a6 => substring(0, 6).padLeft(6);
  String get a15 => substring(0, 15).padLeft(15);
  String get a79 => substring(0, 79).padLeft(79);
}

extension IntArrayFormatExtension on Array<int> {
  String i4(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i4).join(separator);
  String i5(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i5).join(separator);
  String i8(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i8).join(separator);
}

extension DoubleArrayFormatExtension on Array<double> {
  String d12_3(int x) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.d12_3).join();
}
