import 'dart:io';

import 'package:lapack/src/nio.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dblat2.dart';

void main() async {
  final inputFile = File(path.join(currentFilePath(), 'dblat2.in'));
  final nin = Nin(inputFile.openRead());
  final nout = NullNout();
  await dblat2(nin, nout, asyncTestDriver);
}
