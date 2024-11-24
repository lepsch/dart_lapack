// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:dart_lapack/lapack.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dchkaa.dart';

void main() async {
  final (testDriver, nout) = TestDriver.create();
  final inputFile = File(path.join(currentFilePath(), '..', 'dtest.in'));
  final nin = Nin(inputFile.openRead());
  await dchkaa(nin, nout, testDriver);
}
