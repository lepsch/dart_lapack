import 'dart:io';

import 'package:lapack/src/nio.dart';
import 'package:path/path.dart' as path;

import '../test_driver.dart';
import '../utils.dart';
import 'dchkee.dart';

void main() async {
  const inputs = [
    'nep.in',
    'sep.in',
    'se2.in',
    'dsg.in',
    'svd.in',
    'ded.in',
    'dgg.in',
    'dgd.in',
    'dsb.in',
    'dbb.in',
    'glm.in',
    'gqr.in',

    // 'dec.in',

    // 'dbal.in', 'dbak.in', 'dgbal.in', //
    // 'dgbak.in',    'gsv.in', 'csd.in', 'lse.in' //
  ];
  for (final input in inputs) {
    final inputFile = File(path.join(currentFilePath(), '..', input));
    final nin = Nin(inputFile.openRead());
    final nout = NullNout();
    await dchkee(nin, nout, dartTestDriver);
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
