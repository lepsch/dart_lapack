import 'dart:async';
import 'dart:convert';

import 'package:async/async.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

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
    return (await readLine()).trim().split(RegExp(r'\s+'));
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

  Future<void> readArray<T>(Array<T> a, int n) async {
    final parts = await readList();
    if (parts.length < n) throw EOF();
    for (var i = 1; i <= n; i++) {
      a[i] = _parse<T>(parts[i - 1]);
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

    final a = Array<T>(n);
    await readArray(a, n);
    b1.value = a[1];
    if (n == 1) return;
    b2!.value = a[2];
    if (n == 2) return;
    b3!.value = a[3];
    if (n == 3) return;
    b4!.value = a[4];
    if (n == 4) return;
    b5!.value = a[5];
    if (n == 5) return;
    b6!.value = a[6];
    if (n == 6) return;
    b7!.value = a[7];
    if (n == 7) return;
    b8!.value = a[8];
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

  Future<(T1, T2)> read2<T1, T2>() async {
    final parts = await readList();
    if (parts.length < 2) throw EOF();
    return (
      _parse<T1>(parts[0]),
      _parse<T2>(parts[1]),
    );
  }

  T _parse<T>(String s) {
    return switch (T) {
      int => int.parse(s),
      double =>
        double.parse(s.replaceFirstMapped(RegExp('[Dd]([+-])?'), (match) {
          return 'e${match[1] ?? ''}';
        })),
      bool => s.contains(RegExp('[Tt]')),
      String => s,
      _ => throw UnimplementedError(),
    } as T;
  }
}

class Nout {
  final StreamSink<List<int>> _stream;
  Nout(this._stream);

  void println([String? s]) {
    s ??= '';
    _stream.add(utf8.encode('$s\n'));
  }

  Future<void> close() => _stream.close();
}
