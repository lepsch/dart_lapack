// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:lapack/src/nio.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dchkdmd.dart';

void main() async {
  final inputFile = File(path.join(currentFilePath(), '..', 'ddmd.in'));
  final nin = Nin(inputFile.openRead());
  final nout = NullNout();
  await dchkdmd(nin, nout, asyncTestDriver);
}
