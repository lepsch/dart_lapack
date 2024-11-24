// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:async';
import 'dart:convert';

import 'package:async/async.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

class EOF extends Error {}

class Nin {
  final StreamQueue<String> _lineStream;

  Nin(Stream<List<int>> stream)
      : _lineStream = StreamQueue(
            stream.transform(utf8.decoder).transform(const LineSplitter()));

  Future<void> close() async {}

  Future<String> readLine() async {
    try {
      return await _lineStream.next;
    } on StateError catch (_) {
      throw EOF();
    }
  }

  Future<List<String>> readList() async {
    final line = (await readLine()).trim();
    if (line.isEmpty) return [];
    return line.split(RegExp(r'\s+'));
  }

  Future<String> readString() async {
    final line = (await readLine()).trim();
    if (line.startsWith(RegExp('[\'"]'))) {
      final quote = line[0];
      return line
          .substring(1, line.indexOf(RegExp('(?<!$quote)$quote(?!$quote)'), 1))
          .replaceAll('$quote$quote', quote);
    }
    return (await readList()).first;
  }

  Future<List<T>> readData<T>(int n) async {
    final result = <T>[];
    var i = 1, left = n;
    while (left > 0) {
      var parts = await readList();
      for (final part in parts) {
        result.add(_parse<T>(part));
        if (i == n) return result;
        i++;
      }
      left = left - parts.length;
    }
    return result;
  }

  Future<void> readArray<T>(Array<T> a, int n) async {
    final data = await readData<T>(n);
    a.toData().setAll(0, data);
  }

  Future<void> readMatrix<T>(Matrix<T> a, int m, int n) async {
    a = a.having(offset: (x: 0, y: 0));
    for (var i = 0; i < m; i++) {
      final dataIter = (await readData<T>(n)).iterator;
      for (var j = 0; j < n; j++) {
        if (!dataIter.moveNext()) break;
        a[i][j] = dataIter.current;
      }
    }
  }

  Future<int> readInt() async {
    final a = await readData<int>(1);
    return a[0];
  }

  Future<(int, int)> readInt2() async {
    final a = await readData<int>(2);
    return (a[0], a[1]);
  }

  Future<(int, int, int)> readInt3() async {
    final a = await readData<int>(3);
    return (a[0], a[1], a[2]);
  }

  Future<(int, int, int, int)> readInt4() async {
    final a = await readData<int>(4);
    return (a[0], a[1], a[2], a[3]);
  }

  Future<void> readBoxes<T>(
    Box<T> b1, [
    Box<T>? b2,
    Box<T>? b3,
    Box<T>? b4,
    Box<T>? b5,
    Box<T>? b6,
    Box<T>? b7,
    Box<T>? b8,
  ]) async {
    final n = 1 + [b2, b3, b4, b5, b6, b7, b8].nonNulls.length;

    final a = await readData<T>(n);
    b1.value = a[0];
    if (n == 1) return;
    b2!.value = a[1];
    if (n == 2) return;
    b3!.value = a[2];
    if (n == 3) return;
    b4!.value = a[3];
    if (n == 4) return;
    b5!.value = a[4];
    if (n == 5) return;
    b6!.value = a[5];
    if (n == 6) return;
    b7!.value = a[6];
    if (n == 7) return;
    b8!.value = a[7];
  }

  Future<double> readDouble() async {
    final a = await readData<double>(1);
    return a[0];
  }

  Future<(double, double)> readDouble2() async {
    final a = await readData<double>(2);
    return (a[0], a[1]);
  }

  Future<(double, double, double, double)> readDouble4() async {
    final a = await readData<double>(4);
    return (a[0], a[1], a[2], a[3]);
  }

  Future<bool> readBool() async {
    final a = await readData<bool>(1);
    return a[0];
  }

  Future<(T1, T2)> read2<T1, T2>() async {
    final parts = await readList();
    final minLength = null is T1
        ? 0
        : null is T2
            ? 1
            : 2;
    if (parts.length < minLength) throw EOF();

    switch (parts.length) {
      case >= 2:
        return (
          _nonNullableParse<T1>(parts[0]),
          _nonNullableParse<T2>(parts[1]),
        );
      case == 1:
        return (
          _nonNullableParse<T1>(parts[0]),
          null as T2,
        );
      default:
        return (null as T1, null as T2);
    }
  }

  T _nonNullableParse<T>(String value) {
    if (1 is T) return _parse<int>(value) as T;
    if (1.5 is T) return _parse<double>(value) as T;
    if (true is T) return _parse<bool>(value) as T;
    if ('' is T) return _parse<String>(value) as T;
    if (Complex.zero is T) return _parse<Complex>(value) as T;
    throw UnimplementedError();
  }

  T _parse<T>(String s) {
    return switch (T) {
      const (int) => int.parse(s),
      const (double) =>
        double.parse(s.replaceFirstMapped(RegExp('[Dd]([+-])?'), (match) {
          return 'e${match[1] ?? ''}';
        })),
      const (bool) => s.contains(RegExp('[Tt]')),
      const (String) => s,
      const (Complex) => _parseComplex(s),
      _ => throw UnimplementedError(),
    } as T;
  }

  Complex _parseComplex(String s) {
    // Remove parenthesis
    s = s.substring(1, s.length - 1);
    final [real, imaginary] = s.split(',');
    return Complex(_parse<double>(real), _parse<double>(imaginary));
  }
}

abstract class Nout {
  factory Nout(StreamSink<List<int>> stream) = StreamNout;

  void println([String? s]);

  Future<void> close();
}

class StreamNout implements Nout {
  final StreamSink<List<int>> _stream;

  StreamNout(this._stream);

  @override
  void println([String? s]) {
    s ??= '';
    _stream.add(utf8.encode('$s\n'));
  }

  @override
  Future<void> close() => _stream.close();
}

class NoutDelegator<T extends Nout> implements Nout {
  final T nout;

  const NoutDelegator(this.nout);

  @override
  void println([String? s]) => nout.println(s);

  @override
  Future<void> close() => nout.close();
}

class NullNout implements Nout {
  @override
  void println([String? s]) {}

  @override
  Future<void> close() => Future.value();
}
