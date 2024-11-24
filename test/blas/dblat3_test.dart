// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:lapack/blas.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dblat3.dart';

void main() async {
  final inputFile = File(path.join(currentFilePath(), 'dblat3.in'));
  final nin = Nin(inputFile.openRead());
  final nout = NullNout();
  await dblat3(nin, nout, asyncTestDriver);
}
