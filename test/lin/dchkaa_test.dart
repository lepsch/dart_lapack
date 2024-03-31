import 'dart:io';

import 'package:lapack/lapack.dart';
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
