import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';

extension DoubleMatrixExtension on Matrix<double> {
  void debug(String name, int m, int n) {
    dbg('$name(${m.i6}, ${n.i6}) =');
    for (var i = 1; i <= m; i++) {
      var line = '';
      for (var j = 1; j <= n; j++) {
        line += this[i][j].$;
      }
      dbg(line);
    }
  }
}

extension DoubleArrayExtension on Array<double> {
  void debug(String name, int s) {
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += this[i].$;
    }
    dbg(line);
  }
}

extension IntArrayExtension on Array<int> {
  void debug(String name, int s) {
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += this[i].$;
    }
    dbg(line);
  }
}

extension BoolArrayExtension on Array<bool> {
  void debug(String name, int s) {
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += (this[i] ? 1 : 0).i2;
    }
    dbg(line);
  }
}

extension DoubleExtension on double {
  String get $ {
    double v = this == 0 ? 0 : this;
    if (v.abs() < 1e-14) v = 0;
    // return '${v.d10_3.substring(0, 5)}0${v.d10_3.substring(6)}';
    return '${v.d30_20.substring(0, 25)}0${v.d30_20.substring(26)}';
  }
}

extension IntExtension on int {
  String get $ => i4;
}

void debugd1(String name, double d) {
  final a = Array.fromList([d]);
  a.debug(name, 1);
}

void debugd2(String name, double d1, double d2) {
  final a = Array.fromList([d1, d2]);
  a.debug(name, 2);
}

void debugd3(String name, double d1, double d2, double d3) {
  final a = Array.fromList([d1, d2, d3]);
  a.debug(name, 3);
}

void debugd4(String name, double d1, double d2, double d3, double d4) {
  final a = Array.fromList([d1, d2, d3, d4]);
  a.debug(name, 4);
}

void debugd5(
    String name, double d1, double d2, double d3, double d4, double d5) {
  final a = Array.fromList([d1, d2, d3, d4, d5]);
  a.debug(name, 5);
}

void debugd6(String name, double d1, double d2, double d3, double d4, double d5,
    double d6) {
  final a = Array.fromList([d1, d2, d3, d4, d5, d6]);
  a.debug(name, 6);
}

void debugdi(String name, double d, int i) {
  dbg('$name ${d.$}${i.$}');
}

void debugi1(String name, int d) {
  final a = Array.fromList([d]);
  a.debug(name, 1);
}

void debugi2(String name, int d1, int d2) {
  final a = Array.fromList([d1, d2]);
  a.debug(name, 2);
}

void debugi3(String name, int d1, int d2, int d3) {
  final a = Array.fromList([d1, d2, d3]);
  a.debug(name, 3);
}

void debugi4(String name, int d1, int d2, int d3, int d4) {
  final a = Array.fromList([d1, d2, d3, d4]);
  a.debug(name, 4);
}

void debugi5(String name, int d1, int d2, int d3, int d4, int d5) {
  final a = Array.fromList([d1, d2, d3, d4, d5]);
  a.debug(name, 5);
}

void debugi6(String name, int d1, int d2, int d3, int d4, int d5, int d6) {
  final a = Array.fromList([d1, d2, d3, d4, d5, d6]);
  a.debug(name, 6);
}

void debugi7(
    String name, int d1, int d2, int d3, int d4, int d5, int d6, int d7) {
  final a = Array.fromList([d1, d2, d3, d4, d5, d6, d7]);
  a.debug(name, 7);
}

void debugi8(String name, int d1, int d2, int d3, int d4, int d5, int d6,
    int d7, int d8) {
  final a = Array.fromList([d1, d2, d3, d4, d5, d6, d7, d8]);
  a.debug(name, 8);
}

void dbg(String s) {
  if (isInDebugMode) return;
  print(s);
  // stderr.writeln(s);
}

bool get isInDebugMode {
  bool inDebugMode = false;
  assert(inDebugMode = true);
  return inDebugMode;
}
