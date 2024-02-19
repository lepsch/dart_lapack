import 'dart:math';

import 'package:lapack/src/matrix.dart';

extension IntFormatExtension on int {
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

extension DoubleFormatExtension on double {
  // fixed-point
  String f(int w, [int d = 0]) => toStringAsFixed(d).padLeft(w).w(w);
  String get f4_2 => f(4, 2);
  String get f8_2 => f(8, 2);
  String get f12_2 => f(12, 2);
  String get f12_3 => f(12, 3);

  // double
  String d(int w, [int d = 0]) {
    var [num, exponent] = toStringAsExponential(d - 1).split('e');
    num = num.replaceFirstMapped(RegExp('([+-]?)(\\d).'), (match) {
      return '${match[1]}0.${match[2]}';
    });
    final exp = int.parse(exponent) + (this == 0 ? 0 : 1);
    return '$num${exp.abs().toString().length <= 2 ? 'D' : ''}${exp < 0 ? '-' : '+'}${exp.abs().toString().padLeft(2, '0')}'
        .padLeft(w);
  }

  String get d10_3 => d(10, 3);
  String get d12_3 => d(12, 3);
  String get d12_4 => d(12, 4);
  String get d16_6 => d(16, 6);
  String get d36_8 => d(36, 8);

  // exponential
  String e(int w, [int d = 0]) => this.d(w, d).replaceFirst('D', 'E');
  String get e12_3 => e(12, 3);
  String get e15_8 => e(15, 8);
  String get e16_6 => e(16, 6);

  // general
  String g(int w, [int d = 0]) {
    var n = abs();
    const minFp = 0.1;
    final maxFp = pow(10, d);
    if (n != 0) {
      if (n < minFp || n >= maxFp) return e(w, d);
      var prevFp = minFp;
      var nextFp = prevFp * 10;
      if (n >= 10) d -= 1;
      while (d > 0 &&
          !(n >= prevFp && (n < 10 || n == maxFp ? n < nextFp : n <= nextFp))) {
        prevFp = nextFp;
        nextFp *= 10;
        d -= 1;
      }
    } else {
      d -= 1;
    }

    if (d == 0) {
      return '${truncate().toString()}.'.padLeft(w - 4) + ' ' * 4;
    }
    return f(w - 4, d) + ' ' * 4;
  }

  String get g10_3 => g(10, 3);
  String get g12_3 => g(12, 3);
  String get g11_4 => g(11, 4);
  String get g13_6 => g(13, 6);
  String get g16_6 => g(16, 6);
}

extension StringFormatExtension on String {
  String get a1 => substring(0, 1).padLeft(1);
  String get a3 => substring(0, 3).padLeft(3);
  String get a4 => substring(0, 4).padLeft(4);
  String get a6 => substring(0, 6).padLeft(6);
  String get a10 => substring(0, 10).padLeft(10);
  String get a15 => substring(0, 15).padLeft(15);
  String get a79 => substring(0, 79).padLeft(79);

  String w(int w) => length > w ? '*' * w : this;
}

extension IntArrayFormatExtension on Array<int> {
  String i4(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i4).join(separator);
  String i5(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i5).join(separator);
  String i6(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i6).join(separator);
  String i8(int x, [String separator = '']) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.i8).join(separator);
}

extension DoubleArrayFormatExtension on Array<double> {
  String d12_3(int x) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.d12_3).join();
}
