// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart' as test_pkg;

abstract interface class TestDriver {
  const TestDriver();

  static bool get isAsync =>
      Platform.script.normalizePath().path.contains('/dart_test.');

  static (TestDriver, Nout) create() {
    return TestDriver.isAsync
        ? (asyncTestDriver, NullNout())
        : (syncTestDriver, Nout(stdout));
  }

  void group(String description, dynamic Function() body);
  void setUp(dynamic Function() body);
  void tearDown(dynamic Function() body);
  void test(String description, dynamic Function() body, {Object? skip});
  void fail([String message]);
  bool expect(dynamic actual, dynamic matcher, {String? reason});

  void call(String description, dynamic Function() test, {Object? skip}) {
    this.test(description, test, skip: skip);
  }
}

int _errors = 0;

class _SyncTestDriver extends TestDriver {
  const _SyncTestDriver();

  int get errors => _errors;

  @override
  void group(String description, dynamic Function() body) {
    body();
  }

  @override
  void setUp(dynamic Function() body) {
    body();
  }

  @override
  void tearDown(dynamic Function() body) {
    body();
  }

  @override
  void test(String description, dynamic Function() body, {Object? skip}) {
    if (skip is bool) {
      if (skip) return;
    } else if (skip != null) {
      return;
    }
    body();
  }

  @override
  void fail([String message = '']) {
    exit(1);
  }

  @override
  bool expect(dynamic actual, dynamic matcher, {String? reason}) {
    var matchState = <dynamic, dynamic>{};
    if (test_pkg.wrapMatcher(matcher).matches(actual, matchState)) {
      return true;
    }
    _errors++;
    return false;
  }
}

class _AsyncTestDriver extends TestDriver {
  const _AsyncTestDriver();

  @override
  void group(String description, dynamic Function() body) {
    test_pkg.group(description, body);
  }

  @override
  void setUp(dynamic Function() body) {
    test_pkg.setUp(body);
  }

  @override
  void tearDown(dynamic Function() body) {
    test_pkg.tearDown(body);
  }

  @override
  void test(String description, dynamic Function() body, {Object? skip}) {
    test_pkg.test(description, body, skip: skip);
  }

  @override
  void fail([String message = '']) {
    test_pkg.fail(message);
  }

  @override
  bool expect(dynamic actual, dynamic matcher, {String? reason}) {
    test_pkg.expect(actual, matcher, reason: reason);
    return true;
  }
}

const syncTestDriver = _SyncTestDriver();
const asyncTestDriver = _AsyncTestDriver();
