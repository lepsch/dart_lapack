// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/nio.dart';

void alahdg(final Nout IOUNIT, final String PATH) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  String C2;
  int ITYPE = 0;

  // if (IOUNIT <= 0) return;
  C2 = PATH.substring(0, 3);

  // First line describing matrices in this path

  if (lsamen(3, C2, 'GQR')) {
    ITYPE = 1;
    IOUNIT.println('\n $PATH ${_data[9991]}');
  } else if (lsamen(3, C2, 'GRQ')) {
    ITYPE = 2;
    IOUNIT.println('\n $PATH ${_data[9992]}');
  } else if (lsamen(3, C2, 'LSE')) {
    ITYPE = 3;
    IOUNIT.println('\n $PATH ${_data[9993]}');
  } else if (lsamen(3, C2, 'GLM')) {
    ITYPE = 4;
    IOUNIT.println('\n $PATH ${_data[9994]}');
  } else if (lsamen(3, C2, 'GSV')) {
    ITYPE = 5;
    IOUNIT.println('\n $PATH ${_data[9995]}');
  } else if (lsamen(3, C2, 'CSD')) {
    ITYPE = 6;
    IOUNIT.println('\n $PATH ${_data[9996]}');
  }

  // Matrix types

  IOUNIT.println(' Matrix types: ');

  if (ITYPE == 1) {
    IOUNIT.println('''
   1: ${_data[9950]}
   2: ${_data[9952]}
   3: ${_data[9954]}
   4: ${_data[9955]}
   5: ${_data[9956]}
   6: ${_data[9957]}
   7: ${_data[9961]}
   8: ${_data[9962]}
         ''');
  } else if (ITYPE == 2) {
    IOUNIT.println('''
   1: ${_data[9951]}
   2: ${_data[9953]}
   3: ${_data[9954]}
   4: ${_data[9955]}
   5: ${_data[9956]}
   6: ${_data[9957]}
   7: ${_data[9961]}
   8: ${_data[9962]}
         ''');
  } else if (ITYPE == 3) {
    IOUNIT.println('''
   1: ${_data[9950]}
   2: ${_data[9952]}
   3: ${_data[9954]}
   4: ${_data[9955]}
   5: ${_data[9955]}
   6: ${_data[9955]}
   7: ${_data[9955]}
   8: ${_data[9955]}
         ''');
  } else if (ITYPE == 4) {
    IOUNIT.println('''
   1: ${_data[9951]}
   2: ${_data[9953]}
   3: ${_data[9954]}
   4: ${_data[9955]}
   5: ${_data[9955]}
   6: ${_data[9955]}
   7: ${_data[9955]}
   8: ${_data[9955]}
         ''');
  } else if (ITYPE == 5) {
    IOUNIT.println('''
   1: ${_data[9950]}
   2: ${_data[9952]}
   3: ${_data[9954]}
   4: ${_data[9955]}
   5: ${_data[9956]}
   6: ${_data[9957]}
   7: ${_data[9959]}
   8: ${_data[9960]}
         ''');
  } else if (ITYPE == 6) {
    IOUNIT.println('''
   1: ${_data[9963]}
   2: ${_data[9964]}
   3: ${_data[9965]}
         ''');
  }

  // Tests performed

  IOUNIT.println(' Test ratios: ');

  if (ITYPE == 1) {
    // GQR decomposition of rectangular matrices
    IOUNIT.println('''
   1: ${_data[9930]}
   2: ${_data[9931]}
   3: ${_data[9932]}
   4: ${_data[9933]}
         ''');
  } else if (ITYPE == 2) {
    // GRQ decomposition of rectangular matrices
    IOUNIT.println('''
   1" ${_data[9934]}
   2" ${_data[9935]}
   3" ${_data[9932]}
   4" ${_data[9933]}
         ''');
  } else if (ITYPE == 3) {
    // LSE Problem
    IOUNIT.println('''
   1: ${_data[9937]}
   2: ${_data[9938]}
         ''');
  } else if (ITYPE == 4) {
    // GLM Problem

    IOUNIT.println('   1: ${_data[9939]}');
  } else if (ITYPE == 5) {
    // GSVD
    IOUNIT.println('''
   1: ${_data[9940]}
   2: ${_data[9941]}
   3: ${_data[9942]}
   4: ${_data[9943]}
   5: ${_data[9944]}
         ''');
  } else if (ITYPE == 6) {
    // CSD
    IOUNIT.println('''
   ${_data[9910]}
   1: ${_data[9911]}
   2: ${_data[9912]}
   3: ${_data[9913]}
   4: ${_data[9914]}
   5: ${_data[9915]}
   6: ${_data[9916]}
   7: ${_data[9917]}
   8: ${_data[9918]}
   9: ${_data[9919]}
   ${_data[9920]}
  10: ${_data[9921]}
  11: ${_data[9922]}
  12: ${_data[9923]}
  13: ${_data[9924]}
  14: ${_data[9925]}
  15: ${_data[9926]}
''');
  }
}

const _data = {
  9991: 'GQR factorization of general matrices',
  9992: 'GRQ factorization of general matrices',
  9993: 'LSE Problem',
  9994: 'GLM Problem',
  9995: 'Generalized Singular Value Decomposition',
  9996: 'CS Decomposition',

  9950: 'A-diagonal matrix  B-upper triangular',
  9951: 'A-diagonal matrix  B-lower triangular',
  9952: 'A-upper triangular B-upper triangular',
  9953: 'A-lower triangular B-diagonal triangular',
  9954: 'A-lower triangular B-upper triangular',

  9955: 'Random matrices cond(A)=100, cond(B)=10,',

  9956: 'Random matrices cond(A)= sqrt( 0.1/EPS ) cond(B)= sqrt( 0.1/EPS )',
  9957: 'Random matrices cond(A)= 0.1/EPS cond(B)= 0.1/EPS',
  9959: 'Random matrices cond(A)= sqrt( 0.1/EPS ) cond(B)=  0.1/EPS ',
  9960: 'Random matrices cond(A)= 0.1/EPS cond(B)=  sqrt( 0.1/EPS )',

  9961: 'Matrix scaled near underflow limit',
  9962: 'Matrix scaled near overflow limit',
  9963: 'Random orthogonal matrix (Haar measure)',
  9964:
      'Nearly orthogonal matrix with uniformly distributed angles atan2( S, C ) in CS decomposition',
  9965:
      'Random orthogonal matrix with clustered angles atan2( S, C ) in CS decomposition',

  // GQR test ratio

  9930: 'norm( R - Q' ' * A ) / ( min( N, M )*norm( A )* EPS )',
  9931: 'norm( T * Z - Q' ' * B )  / ( min(P,N)*norm(B)* EPS )',
  9932: 'norm( I - Q' '*Q )   / ( N * EPS )',
  9933: 'norm( I - Z' '*Z )   / ( P * EPS )',

  // GRQ test ratio

  9934: 'norm( R - A * Q' ' ) / ( min( N,M )*norm(A) * EPS )',
  9935: 'norm( T * Q - Z' ' * B )  / ( min( P,N ) * norm(B)*EPS )',

  // LSE test ratio

  9937: 'norm( A*x - c )  / ( norm(A)*norm(x) * EPS )',
  9938: 'norm( B*x - d )  / ( norm(B)*norm(x) * EPS )',

  // GLM test ratio

  9939: 'norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*(norm(x)+norm(y))*EPS )',

  // GSVD test ratio

  9940: 'norm( U' ' * A * Q - D1 * R ) / ( min( M, N )*norm( A ) * EPS )',
  9941: 'norm( V' ' * B * Q - D2 * R ) / ( min( P, N )*norm( B ) * EPS )',
  9942: 'norm( I - U' '*U )   / ( M * EPS )',
  9943: 'norm( I - V' '*V )   / ( P * EPS )',
  9944: 'norm( I - Q' '*Q )   / ( N * EPS )',

  // CSD test ratio

  9910: '   2-by-2 CSD',
  9911: 'norm( U1'
      ' * X11 * V1 - C ) / ( max(  P,  Q) * max(norm(I-X'
      '*X),EPS) )',
  9912: 'norm( U1'
      ' * X12 * V2-(-S)) / ( max(  P,M-Q) * max(norm(I-X'
      '*X),EPS) )',
  9913: 'norm( U2'
      ' * X21 * V1 - S ) / ( max(M-P,  Q) * max(norm(I-X'
      '*X),EPS) )',
  9914: 'norm( U2'
      ' * X22 * V2 - C ) / ( max(M-P,M-Q) * max(norm(I-X'
      '*X),EPS) )',
  9915: 'norm( I - U1' '*U1 ) / (   P   * EPS )',
  9916: 'norm( I - U2' '*U2 ) / ( (M-P) * EPS )',
  9917: 'norm( I - V1' '*V1 ) / (   Q   * EPS )',
  9918: 'norm( I - V2' '*V2 ) / ( (M-Q) * EPS )',
  9919: 'principal angle ordering ( 0 or ULP )',
  9920: '   2-by-1 CSD',
  9921: 'norm( U1'
      ' * X11 * V1 - C ) / ( max(  P,  Q) * max(norm(I-X'
      '*X),EPS) )',
  9922: 'norm( U2'
      ' * X21 * V1 - S ) / ( max(  M-P,Q) * max(norm(I-X'
      '*X),EPS) )',
  9923: 'norm( I - U1' '*U1 ) / (   P   * EPS )',
  9924: 'norm( I - U2' '*U2 ) / ( (M-P) * EPS )',
  9925: 'norm( I - V1' '*V1 ) / (   Q   * EPS )',
  9926: 'principal angle ordering ( 0 or ULP )',
};
