// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

// import 'dart:math';

import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/epsilon.dart';
import 'package:lapack/src/intrinsics/maxexponent.dart';
import 'package:lapack/src/intrinsics/minexponent.dart';
import 'package:lapack/src/intrinsics/radix.dart';

// -- LAPACK auxiliary module --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// Standard constants for
// const sp = 1;

// const szero = 0.0;
// const shalf = 0.5;
// const sone = 1.0;
// const stwo = 2.0;
// const sthree = 3.0;
// const sfour = 4.0;
// const seight = 8.0;
// const sten = 10.0;
// const czero = Complex.zero;
// const chalf = Complex(0.5, 0.0);
// const cone = Complex.one;
// const sprefix = 'S';
// const cprefix = 'C';

// Scaling constants
// const sulp = epsilon(0.0);
// const seps = sulp * 0.5;
// const ssafmin =
//     pow(real(radix(0.0), sp), max(minexponent(0.0) - 1, 1 - maxexponent(0.0)));
// const ssafmax = sone / ssafmin;
// const ssmlnum = ssafmin / sulp;
// const sbignum = ssafmax * sulp;
// const srtmin = sqrt(ssmlnum);
// const srtmax = sqrt(sbignum);

// Blue's scaling constants
// const stsml = pow(real(radix(0.0), sp), ceiling((minexponent(0.0) - 1) * 0.5));
// const stbig = pow(
//     real(radix(0.0), sp), floor((maxexponent(0.0) - digits(0.0) + 1) * 0.5));
// // ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
// // The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
// const sssml =
//     pow(real(radix(0.0), sp), (-floor((minexponent(0.0) - digits(0.0)) * 0.5)));
// // sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
// const ssbig = pow(real(radix(0.0), sp),
//     (-ceiling((maxexponent(0.0) + digits(0.0) - 1) * 0.5)));

// Standard constants for
const dp = 1;

const dzero = 0.0;
const dhalf = 0.5;
const done = 1.0;
const dtwo = 2.0;
const dthree = 3.0;
const dfour = 4.0;
const deight = 8.0;
const dten = 10.0;
const zzero = Complex.zero;
const zhalf = Complex(0.5, 0.0);
const zone = Complex.one;
const dprefix = 'D';
const zprefix = 'Z';

// Scaling constants
final dulp = epsilon(0.0);
final deps = dulp * 0.5;
final dsafmin =
    pow(radix(0.0), max(minexponent(0.0) - 1, 1 - maxexponent(0.0))).toDouble();
final dsafmax = done / dsafmin;
final dsmlnum = dsafmin / dulp;
final dbignum = dsafmax * dulp;
final drtmin = sqrt(dsmlnum);
final drtmax = sqrt(dbignum);

// Blue's scaling constants
final dtsml = pow(radix(0.0), ((minexponent(0.0) - 1) * 0.5).ceil());
final dtbig =
    pow(radix(0.0), ((maxexponent(0.0) - digits(0.0) + 1) * 0.5).floor());
// ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
// The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
final dssml =
    pow(radix(0.0), -((minexponent(0.0) - digits(0.0)) * 0.5).floor());
// sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
final dsbig =
    pow(radix(0.0), -((maxexponent(0.0) + digits(0.0) - 1) * 0.5).ceil());
