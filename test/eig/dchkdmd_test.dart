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
