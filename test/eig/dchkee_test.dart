import 'dart:io';

import 'package:lapack/src/nio.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dchkee.dart';

void main() async {
  final (testDriver, nout) = TestDriver.isTesting
      ? (asyncTestDriver, NullNout())
      : (syncTestDriver, Nout(stdout));

  const inputs = [
    'csd.in',
    'dbak.in',
    'dbal.in',
    'dbb.in',
    'dec.in',
    'ded.in',
    'dgbak.in',
    'dgbal.in',
    'dgd.in',
    'dgg.in',
    'dsb.in',
    'dsg.in',
    'glm.in',
    'gqr.in',
    'gsv.in',
    'lse.in',
    'nep.in',
    'se2.in',
    'sep.in',
    'svd.in',
  ];
  for (final input in inputs) {
    final inputFile = File(path.join(currentFilePath(), '..', input));
    final nin = Nin(inputFile.openRead());
    await dchkee(nin, nout, testDriver);
  }

// dchkdmd.dart
// 'ddmd.in',

// Z
// 'nep.in',
// 'sep.in',
// 'se2.in',
// 'svd.in',
// 'zec.in',
// 'zed.in',
// 'zgg.in',
// 'zgd.in',
// 'zsb.in',
// 'zsg.in',
// 'zbal.in',
// 'zbak.in',
// 'zgbal.in',
// 'zgbak.in',
// 'zbb.in',
// 'glm.in',
// 'gqr.in',
// 'gsv.in',
// 'csd.in',
// 'lse.in',
// 'zdmd.in',
}
