import 'dart:io';

import 'package:test/test.dart' as test_pkg;

abstract interface class TestDriver {
  const TestDriver();

  void group(String description, dynamic Function() body);
  void test(String description, dynamic Function() body, {Object? skip});
  void fail([String message]);
  bool expect(dynamic actual, dynamic matcher);

  void call(String description, dynamic Function() test, {Object? skip}) {
    this.test(description, test);
  }
}

int _errors = 0;

class _LapackTestDriver extends TestDriver {
  const _LapackTestDriver();

  int get errors => _errors;

  @override
  void group(String description, dynamic Function() body) {
    body();
  }

  @override
  void test(String description, dynamic Function() body, {Object? skip}) {
    if (skip != null) return;
    body();
  }

  @override
  void fail([String message = '']) {
    exit(1);
  }

  @override
  bool expect(dynamic actual, dynamic matcher) {
    var matchState = <dynamic, dynamic>{};
    if (test_pkg.wrapMatcher(matcher).matches(actual, matchState)) {
      return true;
    }
    _errors++;
    return false;
  }
}

class _DartTestDriver extends TestDriver {
  const _DartTestDriver();

  @override
  void group(String description, dynamic Function() body) {
    test_pkg.group(description, body);
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
  bool expect(dynamic actual, dynamic matcher) {
    test_pkg.expect(actual, matcher);
    return true;
  }
}

const lapackTestDriver = _LapackTestDriver();
const dartTestDriver = _DartTestDriver();