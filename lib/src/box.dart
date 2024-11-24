// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

class Box<T> {
  T value;
  Box(this.value);

  @override
  String toString() => value.toString();
}

class DelegatingBox<T> implements Box<T> {
  final T Function() _get;
  final void Function(T) _set;

  const DelegatingBox(this._get, this._set);

  @override
  T get value => _get();
  @override
  set value(T value) => _set(value);

  @override
  String toString() => value.toString();
}
