import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';

abstract interface class Array<T> implements Box<T> {
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
  int maxloc(int start, int end, {int dim = 1});

  Array<T> dim([int? ld]);
}

class Matrix<T> implements Box<T>  {
  final Array<T> _entries;
  final int _ld;

  Matrix(int m, int n)
      : _ld = m,
        _entries = _Array<T>(m * n, ld: m);

  Matrix.fromList(List<List<T>> list)
      : _entries = _Array<T>.fromList(
          [
            for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
              for (var i = 0; i < list.length; i++) list[i][j],
            ],
          ],
          ld: list.length,
        ),
        _ld = list.length;

  Matrix.fromSlice(this._entries, this._ld);

  Matrix<T> call(int i, int j, {int? ld}) {
    var entries = _entries(j, offset: i - 1);
    if (ld != null) {
      entries = _Array.fromSlice(entries.toRawList(), ld: ld);
    }
    return Matrix.fromSlice(entries, ld ?? _ld);
  }

  Array<T> operator [](int i) {
    return _MatrixArrayAdapter(this, i);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j) => this[i].box(j);

  Array<T> asArray() => _Array.fromSlice(_entries.toRawList());

  Matrix<T> dim(int ld) => this(1, 1, ld: ld);

  @override
  T get value => this[1][1];

  @override
  set value(T value) => this[1][1] = value;
}

class _MatrixArrayAdapter<T> implements Array<T> {
  final int i;
  final Array<T> _entries;
  final Matrix<T> _m;

  _MatrixArrayAdapter(this._m, this.i) : _entries = _m._entries;

  @override
  Array<T> slice(int j, {int offset = 0}) {
    return _entries.slice(j, offset: i - 1 + offset);
  }

  @override
  Array<T> call(int j, {int offset = 0}) {
    return _entries(j, offset: i - 1 + offset);
  }

  @override
  Box<T> box(int j) {
    return _entries(j, offset: i - 1).box(1);
  }

  @override
  List<T> toRawList() {
    return _entries(1, offset: i - 1).toRawList();
  }

  @override
  Matrix<T> asMatrix(int ld) {
    return _entries(1, offset: i - 1).asMatrix(ld);
  }

  @override
  T maxval(int start, int end) {
    return _entries(1, offset: i - 1).maxval(start, end);
  }

  @override
  int maxloc(int start, int end, {int dim = 1}) {
    return _entries(1, offset: i - 1).maxloc(start, end, dim: dim);
  }

  @override
  T operator [](int j) {
    return _entries(j, offset: i - 1)[1];
  }

  @override
  void operator []=(int j, T value) {
    _entries(j, offset: i - 1)[1] = value;
  }

  @override
  Array<T> dim([int? ld]) => this;

  @override
  T get value => this[1];

  @override
  set value(T value) => this[1] = value;
}

class _Array<T> implements Array<T> {
  final int offset;
  final int _ld;
  final List<T> _elements;

  _Array(int length, {this.offset = 0, int ld = 1})
      : _ld = ld,
        _elements = switch (T) {
          double => Float64List(length),
          int => Int64List(length),
          bool => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromList(List<T> list, {this.offset = 0, int ld = 1})
      : _ld = ld,
        _elements = switch (T) {
          double => Float64List.fromList(list as List<double>),
          int => Int64List.fromList(list as List<int>),
          bool => [...list],
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromSlice(this._elements, {this.offset = 0, int ld = 1}) : _ld = ld;

  @override
  Array<T> slice(int index, {int offset = 0}) {
    return _Array.fromSlice(
      switch (T) {
        double => (_elements as Float64List).buffer.asFloat64List(
              (_elements as Float64List).offsetInBytes +
                  ((index - 1) * _ld + offset) *
                      (_elements as Float64List).elementSizeInBytes,
            ),
        int => (_elements as Int64List).buffer.asInt64List(
              (_elements as Int64List).offsetInBytes +
                  ((index - 1) * _ld + offset) *
                      (_elements as Int64List).elementSizeInBytes,
            ),
        bool => _elements is ListSlice
            ? _elements.slice((index - 1) * _ld + offset, _elements.length)
            : ListSlice(
                _elements, (index - 1) * _ld + offset, _elements.length),
        _ => throw UnimplementedError(),
      } as List<T>,
      ld: _ld,
    );
  }

  @override
  Array<T> call(int index, {int offset = 0}) {
    return slice(index, offset: offset);
  }

  @override
  T operator [](int index) {
    return _elements[(index - 1) * _ld + offset];
  }

  @override
  void operator []=(int index, T value) {
    _elements[(index - 1) * _ld + offset] = value;
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
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) value = this[i] as num;
        }
        return value as T;

      case bool:
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return true as T;
        }
        return false as T;
      default:
        throw UnimplementedError();
    }
  }

  @override
  int maxloc(int start, int end, {int dim = 1}) {
    switch (T) {
      case double:
      case int:
        int loc = 0;
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) {
            value = this[i] as num;
            loc = i;
          }
        }
        return loc;

      case bool:
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return i;
        }
        return 0;
      default:
        throw UnimplementedError();
    }
  }

  @override
  Array<T> dim([int? ld]) => this;

  @override
  T get value => _elements[offset];

  @override
  set value(T value) => _elements[offset] = value;
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
  final (int, int) _ld;

  Matrix3d(int m, int n, int o)
      : _ld = (m, n),
        _entries = _Array<T>(m * n * o, ld: m);

  Matrix3d.fromList(List<List<List<T>>> list)
      : _entries = _Array<T>.fromList(
          [
            for (var k = 0; k < (list.firstOrNull ?? []).length; k++) ...[
              for (var i = 0; i < (list.firstOrNull ?? []).length; i++) ...[
                for (var j = 0; j < list.length; j++) list[i][j][k],
              ],
            ],
          ],
          ld: list.length,
        ),
        _ld = (list.length, list[0].length);

  Matrix3d.fromSlice(this._entries, this._ld);

  Matrix3d<T> call(int i, int j, int k) {
    final ld1 = _ld.$1 - i + 1;
    var entries = _entries(1 + (i - 1) * ld1, offset: j - 1)(k);
    return Matrix3d.fromSlice(entries, _ld);
  }

  Matrix<T> operator [](int i) {
    return _entries(1, offset: (i - 1) * _ld.$1).asMatrix(_ld.$1 * _ld.$2);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j, int k) => this[i][j].box(k);

  Array<T> asArray() => _Array.fromSlice(_entries.toRawList());
}
