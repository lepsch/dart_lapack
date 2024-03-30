import 'dart:io';

import 'package:lapack/blas.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'zblat3.dart';

void main() async {
  final inputFile = File(path.join(currentFilePath(), 'zblat3.in'));
  final nin = Nin(inputFile.openRead());
  final nout = NullNout();
  await zblat3(nin, nout, asyncTestDriver);
}
