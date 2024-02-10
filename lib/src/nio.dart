import 'dart:async';
import 'dart:convert';

import 'package:async/async.dart';
import 'package:lapack/src/matrix.dart';

class EOF extends Error {}

class Nin {
  final StreamQueue<String> _lineStream;

  Nin(Stream<List<int>> stream)
      : _lineStream = StreamQueue(
            stream.transform(utf8.decoder).transform(const LineSplitter()));

  Future<String> readLine() async {
    return await _lineStream.next;
  }

  Future<void> readArray<T>(Array<T> a, int n) async {
    List<String> parts;
    try {
      parts = (await readLine()).split(RegExp(r'\s+'));
    } on StateError catch (_) {
      throw EOF();
    }
    if (parts.length < n) throw EOF();
    for (var i = 1; i <= n; i++) {
      a[i] = switch (T) {
        int => int.parse(parts[i - 1]),
        double => double.parse(parts[i - 1]),
        bool => parts[i - 1].contains(RegExp('Tt')),
        _ => throw UnimplementedError(),
      } as T;
    }
  }

  Future<void> readMatrix<T>(Matrix<T> a, int m, int n) async {
    for (var i = 1; i <= m; i++) {
      await readArray<T>(a[i], n);
    }
  }

  Future<int> readInt() async {
    final a = Array<int>(1);
    await readArray(a, 1);
    return a[1];
  }

  Future<(int, int)> readInt2() async {
    final a = Array<int>(2);
    await readArray(a, 2);
    return (a[1], a[2]);
  }

  Future<(int, int, int)> readInt3() async {
    final a = Array<int>(3);
    await readArray(a, 3);
    return (a[1], a[2], a[3]);
  }

  Future<(int, int, int, int)> readInt4() async {
    final a = Array<int>(4);
    await readArray(a, 4);
    return (a[1], a[2], a[3], a[4]);
  }

  Future<double> readDouble() async {
    final a = Array<double>(1);
    await readArray(a, 1);
    return a[1];
  }

  Future<(double, double)> readDouble2() async {
    final a = Array<double>(2);
    await readArray(a, 2);
    return (a[1], a[2]);
  }

  Future<(double, double, double, double)> readDouble4() async {
    final a = Array<double>(4);
    await readArray(a, 4);
    return (a[1], a[2], a[3], a[4]);
  }

  Future<bool> readBool() async {
    final a = Array<bool>(1);
    await readArray(a, 1);
    return a[1];
  }
}

class Nout {
  StreamSink<List<int>> _stream;
  Nout(this._stream);

  void println([String? s]) {
    s ??= '';
    _stream.add(utf8.encode('$s\n'));
  }
}
