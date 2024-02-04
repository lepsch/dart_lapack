import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';

abstract interface class Array<T> {
  factory Array(int length, {int offset = 0}) {
    return _Array<T>(length, offset: offset);
  }

  factory Array.fromList(List<T> list, {int offset = 0}) {
    return _Array.fromList(list, offset: offset);
  }

  factory Array.fromSlice(List<T> elements, {int offset = 0}) {
    return _Array.fromSlice(elements, offset: offset);
  }

  Array<T> slice(int index, {int offset = 0});

  Array<T> call(int index, {int offset = 0});

  T operator [](int index);
  void operator []=(int index, T value);

  Box<T> box(index);

  List<T> toRawList();
}

class Matrix<T> {
  final Array<T> _entries;
  final ({int m, int n}) dimension;
  final ({int x, int y}) offset;

  Matrix({required int m, required int n, this.offset = (x: 0, y: 0)})
      : dimension = (m: m, n: n),
        _entries = _Array<T>(m * n, offset: offset.x, ld: m);

  Matrix.fromList(List<List<T>> list, {this.offset = (x: 0, y: 0)})
      : _entries = _Array<T>.fromList([
          for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
            for (var i = 0; i < list.length; i++) list[i][j]
          ]
        ], ld: list.length),
        dimension = (m: list.length, n: list.isEmpty ? 0 : list[0].length);

  Matrix.fromSlice(this._entries, this.dimension, {this.offset = (x: 0, y: 0)});

  Matrix<T> call(int i, int j, [({int x, int y}) offset = (x: 0, y: 0)]) {
    return Matrix.fromSlice(_entries, dimension,
        offset: (x: offset.x + j - 1, y: offset.y + i - 1));
  }

  Array<T> operator [](int i) {
    return _entries(1 + offset.x, offset: offset.y + i - 1);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j) => this[i].box(j);
}

class _Array<T> implements Array<T> {
  final int offset;
  final int ld;
  final List<T> _elements;

  _Array(int length, {this.offset = 0, this.ld = 1})
      : _elements = switch (T) {
          double => Float64List(length),
          int => Int64List(length),
          bool => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromList(List<T> list, {this.offset = 0, this.ld = 1})
      : _elements = switch (T) {
          double => Float64List.fromList(list as List<double>),
          int => Int64List.fromList(list as List<int>),
          bool => [...list],
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromSlice(this._elements, {this.offset = 0, this.ld = 1});

  @override
  Array<T> slice(int index, {int offset = 0}) {
    return _Array.fromSlice(
      switch (T) {
        double => (_elements as Float64List).buffer.asFloat64List(
            (_elements as Float64List).offsetInBytes +
                ((index - 1) * ld + offset) *
                    (_elements as Float64List).elementSizeInBytes),
        int => (_elements as Int64List).buffer.asInt64List(
            (_elements as Int64List).offsetInBytes +
                ((index - 1) * ld + offset) *
                    (_elements as Int64List).elementSizeInBytes),
        bool => _elements is ListSlice
            ? _elements.slice((index - 1) * ld + offset, _elements.length)
            : ListSlice(_elements, (index - 1) * ld + offset, _elements.length),
        _ => throw UnimplementedError(),
      } as List<T>,
      ld: ld,
    );
  }

  @override
  Array<T> call(int index, {int offset = 0}) {
    return slice(index, offset: offset);
  }

  @override
  T operator [](int index) {
    return _elements[(index - 1) * ld + offset];
  }

  @override
  void operator []=(int index, T value) {
    _elements[(index - 1) * ld + offset] = value;
  }

  @override
  Box<T> box(index) {
    return _ArrayElementBox(this, index);
  }

  @override
  List<T> toRawList() {
    return _elements;
  }
}

class _ArrayElementBox<T> implements Box<T> {
  final Array<T> _array;
  final int _index;
  _ArrayElementBox(this._array, this._index);

  @override
  T get value => _array[_index];
  @override
  set value(T value) => _array[_index] = value;

  @override
  String toString() => value.toString();
}
