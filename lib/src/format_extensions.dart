import 'dart:math';

import 'package:lapack/src/matrix.dart';

extension IntFormatExtension on int {
  String i(int w) => toString().padLeft(w).w(w);

  String get i1 => i(1);
  String get i2 => i(2);
  String get i3 => i(3);
  String get i4 => i(4);
  String get i5 => i(5);
  String get i6 => i(6);
  String get i7 => i(7);
  String get i8 => i(8);
  String get i12 => i(12);
  String get i15 => i(15);
  String get i36 => i(36);
}

extension DoubleFormatExtension on double {
  // fixed-point
  String f(int w, [int d = 0]) => toStringAsFixed(d).padLeft(w).w(w);
  String get f4_1 => f(4, 1);
  String get f4_2 => f(4, 2);
  String get f5_0 => f(5, 0);
  String get f6_1 => f(6, 1);
  String get f8_2 => f(8, 2);
  String get f12_2 => f(12, 2);
  String get f12_3 => f(12, 3);

  // double
  String d(int w, [int d = 0]) {
    if (isNaN) return 'NaN'.padLeft(w);

    var [num, exponent] = toStringAsExponential(d - 1).split('e');
    num = num.replaceFirstMapped(RegExp('([+-]?)(\\d).'), (match) {
      return '${match[1]}0.${match[2]}';
    });
    final exp = int.parse(exponent) + (this == 0 ? 0 : 1);
    return '$num${exp.abs().toString().length <= 2 ? 'D' : ''}${exp < 0 ? '-' : '+'}${exp.abs().toString().padLeft(2, '0')}'
        .padLeft(w);
  }

  String get d9_1 => d(9, 1);
  String get d10_3 => d(10, 3);
  String get d12_3 => d(12, 3);
  String get d12_4 => d(12, 4);
  String get d16_6 => d(16, 6);
  String get d20_6 => d(20, 6);
  String get d30_20 => d(30, 20);
  String get d36_8 => d(36, 8);

  // exponential
  String e(int w, [int d = 0, int e = 0]) =>
      this.d(w, d).replaceFirst('D', 'E');
  String get e8_1 => e(8, 1);
  String get e12_3 => e(12, 3);
  String get e15_8 => e(15, 8);
  String get e16_6 => e(16, 6);
  String get e24_16 => e(24, 16);
  String get e24_16e3 => e(24, 16, 3);

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
  String get g18_6 => g(18, 6);

  String get p => sp;
  String get sp => this >= 0 ? '+' : '-';
  String get ss => this >= 0 ? '' : '-';
}

extension StringFormatExtension on String {
  String a(int w) => padLeft(w).substring(0, w);

  String get a1 => a(1);
  String get a2 => a(2);
  String get a3 => a(3);
  String get a4 => a(4);
  String get a6 => a(6);
  String get a7 => a(7);
  String get a10 => a(10);
  String get a15 => a(15);
  String get a17 => a(17);
  String get a79 => a(79);

  String w(int w) => length > w ? '*' * w : this;
}

extension BoolFormatExtension on bool {
  String get l1 => (this ? 'T' : 'F');
  String get l2 => (this ? ' T' : ' F');
}

extension StringIterableFormatExtension on Iterable<String> {
  String a(int w, [int? len, String separator = '']) =>
      (len != null ? take(len) : this).map((s) => s.a(w)).join(separator);

  String a1([int? len, String separator = '']) => a(1, len, separator);
}

extension IntIterableFormatExtension on Iterable<int> {
  String i(int w, [int? len, String separator = '']) =>
      (len != null ? take(len) : this).map((n) => n.i(w)).join(separator);

  String i3([int? len, String separator = '']) => i(3, len, separator);
  String i5([int? len, String separator = '']) => i(5, len, separator);
  String i36([int? len, String separator = '']) => i(36, len, separator);
}

extension DoubleIterableFormatExtension on Iterable<double> {
  String d(int w, [int d = 0,int? len, String separator = '']) =>
      (len != null ? take(len) : this).map((n) => n.d(w)).join(separator);

  String d3([int? len, String separator = '']) => d(3, 0, len, separator);
  String d5([int? len, String separator = '']) => d(5, 0, len, separator);
  String d12_4([int? len, String separator = '']) => d(12,4, len, separator);
  String d36_8([int? len, String separator = '']) => d(36, 8, len, separator);
}

extension IntArrayFormatExtension on Array<int> {
  String i(int w, int len, [String separator = '']) =>
      [for (var i = 1; i <= len; i++) this[i]].i(w, len);

  String i4(int len, [String separator = '']) => i(4, len, separator);
  String i5(int len, [String separator = '']) => i(5, len, separator);
  String i6(int len, [String separator = '']) => i(6, len, separator);
  String i8(int len, [String separator = '']) => i(8, len, separator);
}

extension DoubleArrayFormatExtension on Array<double> {
  String d12_3(int x) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.d12_3).join();
  String f6_1(int x) =>
      [for (var i = 1; i <= x; i++) this[i]].map((n) => n.f6_1).join();
}
