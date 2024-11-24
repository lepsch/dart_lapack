// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

// import 'dart:io';

import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
// ignore: depend_on_referenced_packages
import 'package:stack_trace/stack_trace.dart';

var _debug = true;

void debugOff() {
  _debug = false;
}

void debugOn() {
  _debug = true;
}

extension DoubleMatrixExtension on Matrix<double> {
  void debug(String name, int m, int n) {
    if (!_debug) return;
    _print('$name(${m.i6}, ${n.i6}) =');
    for (var i = 1; i <= m; i++) {
      var line = '';
      for (var j = 1; j <= n; j++) {
        line += this[i][j].$;
      }
      _print(line);
    }
  }
}

extension DoubleArrayExtension on Array<double> {
  void debug(String name, int s) {
    if (!_debug) return;
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += this[i].$;
    }
    _print(line);
  }
}

extension IntArrayExtension on Array<int> {
  void debug(String name, int s) {
    if (!_debug) return;
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += this[i].$;
    }
    _print(line);
  }
}

extension BoolArrayExtension on Array<bool> {
  void debug(String name, int s) {
    if (!_debug) return;
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += (this[i] ? 1 : 0).i2;
    }
    _print(line);
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
  if (!_debug) return;
  final a = Array.fromList([d]);
  a.debug(name, 1);
}

void debugd2(String name, double d1, double d2) {
  if (!_debug) return;
  final a = Array.fromList([d1, d2]);
  a.debug(name, 2);
}

void debugd3(String name, double d1, double d2, double d3) {
  if (!_debug) return;
  final a = Array.fromList([d1, d2, d3]);
  a.debug(name, 3);
}

void debugd4(String name, double d1, double d2, double d3, double d4) {
  if (!_debug) return;
  final a = Array.fromList([d1, d2, d3, d4]);
  a.debug(name, 4);
}

void debugd5(
    String name, double d1, double d2, double d3, double d4, double d5) {
  if (!_debug) return;
  final a = Array.fromList([d1, d2, d3, d4, d5]);
  a.debug(name, 5);
}

void debugd6(String name, double d1, double d2, double d3, double d4, double d5,
    double d6) {
  if (!_debug) return;
  final a = Array.fromList([d1, d2, d3, d4, d5, d6]);
  a.debug(name, 6);
}

void debugdi(String name, double d, int i) {
  if (!_debug) return;
  debug('$name ${d.$}${i.$}');
}

void debugi1(String name, int i) {
  if (!_debug) return;
  final a = Array.fromList([i]);
  a.debug(name, 1);
}

void debugi2(String name, int i1, int i2) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2]);
  a.debug(name, 2);
}

void debugi3(String name, int i1, int i2, int i3) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3]);
  a.debug(name, 3);
}

void debugi4(String name, int i1, int i2, int i3, int i4) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3, i4]);
  a.debug(name, 4);
}

void debugi5(String name, int i1, int i2, int i3, int i4, int i5) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3, i4, i5]);
  a.debug(name, 5);
}

void debugi6(String name, int i1, int i2, int i3, int i4, int i5, int i6) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3, i4, i5, i6]);
  a.debug(name, 6);
}

void debugi7(
    String name, int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3, i4, i5, i6, i7]);
  a.debug(name, 7);
}

void debugi8(String name, int i1, int i2, int i3, int i4, int i5, int i6,
    int i7, int i8) {
  if (!_debug) return;
  final a = Array.fromList([i1, i2, i3, i4, i5, i6, i7, i8]);
  a.debug(name, 8);
}

void debug([String? s]) {
  if (!_debug) return;
  if (s == null) {
    final (:file, :function, :line) = getCallerInfo();
    print('$file:$line - $function');
    return;
  }
  _print(s);
}

void _print(String s) {
  // if (isInDebugMode) return;
  print(s);
  // stderr.writeln(s);
}

bool get isInDebugMode {
  bool inDebugMode = false;
  assert(inDebugMode = true);
  return inDebugMode;
}

({String file, String function, int line}) getCallerInfo([int depth = 1]) {
  final trace = Trace.from(StackTrace.current).terse;
  final prevFrame = trace.frames[depth + 1];
  return (
    line: prevFrame.line ?? 0,
    file: prevFrame.uri.toString(),
    function: prevFrame.member ?? '',
  );
}
