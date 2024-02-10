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

  Box<T> box(int index);

  List<T> toRawList();

  Matrix<T> asMatrix(int ld);

  T maxval(int start, int end);
}

class Matrix<T> {
  final Array<T> _entries;
  final int _ld;
  final ({int x, int y}) offset;

  Matrix(int m, int n, {this.offset = (x: 0, y: 0)})
      : _ld = m,
        _entries = _Array<T>(m * n, offset: offset.x, ld: m);

  Matrix.fromList(List<List<T>> list, {this.offset = (x: 0, y: 0)})
      : _entries = _Array<T>.fromList(
          [
            for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
              for (var i = 0; i < list.length; i++) list[i][j],
            ],
          ],
          ld: list.length,
        ),
        _ld = list.length;

  Matrix.fromSlice(this._entries, this._ld, {this.offset = (x: 0, y: 0)});

  Matrix<T> call(int i, int j, [({int x, int y}) offset = (x: 0, y: 0)]) {
    return Matrix.fromSlice(_entries, _ld,
        offset: (x: offset.x + j - 1, y: offset.y + i - 1));
  }

  Array<T> operator [](int i) {
    return _entries(1 + offset.x, offset: offset.y + i - 1);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j) => this[i].box(j);

  Array<T> asArray() =>
      _Array.fromSlice(_entries.toRawList(), offset: offset.x * _ld + offset.y);
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
                        (_elements as Float64List).elementSizeInBytes,
              ),
          int => (_elements as Int64List).buffer.asInt64List(
                (_elements as Int64List).offsetInBytes +
                    ((index - 1) * ld + offset) *
                        (_elements as Int64List).elementSizeInBytes,
              ),
          bool => _elements is ListSlice
              ? _elements.slice((index - 1) * ld + offset, _elements.length)
              : ListSlice(
                  _elements, (index - 1) * ld + offset, _elements.length),
          _ => throw UnimplementedError(),
        } as List<T>,
        ld: ld);
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
  Box<T> box(int index) {
    return _ArrayElementBox(this, index);
  }

  @override
  List<T> toRawList() {
    return _elements;
  }

  @override
  Matrix<T> asMatrix(int ld) {
    return Matrix.fromSlice(
        _Array.fromSlice(_elements, offset: offset, ld: ld), ld);
  }

  @override
  T maxval(int start, int end) {
    switch (T) {
      case double:
      case int:
        num value = 0;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[start] as num)) value = this[start] as num;
        }
        return value as T;

      case bool:
        for (var i = start + 1; i <= end; i++) {
          if (this[start] as bool) return true as T;
        }
        return false as T;
      default:
        throw UnimplementedError();
    }
    ;
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

class Matrix3d<T> {
  final Array<T> _entries;
  final int _ld;
  final ({int x, int y, int z}) offset;

  Matrix3d(int m, int n, int o, {this.offset = (x: 0, y: 0, z: 0)})
      : _ld = m,
        _entries = _Array<T>(m * n * o, offset: offset.x, ld: m);

  Matrix3d.fromList(List<List<T>> list, {this.offset = (x: 0, y: 0, z: 0)})
      : _entries = _Array<T>.fromList(
          [
            for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
              for (var i = 0; i < list.length; i++) list[i][j],
            ],
          ],
          ld: list.length,
        ),
        _ld = list.length;

  Matrix3d.fromSlice(this._entries, this._ld,
      {this.offset = (x: 0, y: 0, z: 0)});

  Matrix3d<T> call(int i, int j, int k,
      [({int x, int y, int z}) offset = (x: 0, y: 0, z: 0)]) {
    return Matrix3d.fromSlice(_entries, _ld, offset: (
      x: offset.x + j - 1,
      y: offset.y + i - 1,
      z: offset.z + k - 1
    ));
  }

  Matrix<T> operator [](int i) {
    return Matrix(1, 1);
    // _entries(1 + offset.x, offset: offset.y + i - 1);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j, int k) => this[i][j].box(k);

  Array<T> asArray() =>
      _Array.fromSlice(_entries.toRawList(), offset: offset.x * _ld + offset.y);
}
