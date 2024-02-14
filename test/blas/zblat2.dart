import 'common.dart';

      void main() {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 17 ;
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      int                NMAX, INCMAX;
      const              NMAX = 65, INCMAX = 2 ;
      int                NINMAX, NIDMAX, NKBMAX, NALMAX, NBEMAX;
      const              NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      double             EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB, NOUT, NTRA;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR;
      String             TRANS;
      String             SNAMET;
      String             SNAPS, SUMMRY;
      // .. Local Arrays ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BET( NBEMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( 2*NMAX );
      double             G( NMAX );
      int                IDIM( NIDMAX ), INC( NINMAX ), KB( NKBMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      //- double             DDIFF;
      //- bool               LZE;
      // EXTERNAL DDIFF, LZE
      // .. External Subroutines ..
      // EXTERNAL ZCHK1, ZCHK2, ZCHK3, ZCHK4, ZCHK5, ZCHK6, ZCHKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      // int                infoc.INFOT, infoc.NOUTC;
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // .. Common blocks ..
      // COMMON /infoc/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // COMMON /srnamc/srnamc.SRNAMT
      // .. Data statements ..
      const SNAMES = ['ZGEMV ', 'ZGBMV ', 'ZHEMV ', 'ZHBMV ', 'ZHPMV ', 'ZTRMV ', 'ZTBMV ', 'ZTPMV ', 'ZTRSV ', 'ZTBSV ', 'ZTPSV ', 'ZGERC ', 'ZGERU ', 'ZHER  ', 'ZHPR  ', 'ZHER2 ', 'ZHPR2 '];
      // .. Executable Statements ..

      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SUMMRY;
      READ( NIN, FMT = * )NOUT;
      OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' );
      infoc.NOUTC = NOUT;

      // Read name and unit number for snapshot output file and open file.

      READ( NIN, FMT = * )SNAPS;
      READ( NIN, FMT = * )NTRA;
      TRACE = NTRA >= 0;
      if ( TRACE ) {
         OPEN( NTRA, FILE = SNAPS, STATUS = 'UNKNOWN' );
      }
      // Read the flag that directs rewinding of the snapshot file.
      READ( NIN, FMT = * )REWI;
      REWI = REWI && TRACE;
      // Read the flag that directs stopping on any failure.
      READ( NIN, FMT = * )SFATAL;
      // Read the flag that indicates whether error exits are to be tested.
      READ( NIN, FMT = * )TSTERR;
      // Read the threshold value of the test ratio
      READ( NIN, FMT = * )THRESH;

      // Read and check the parameter values for the tests.

      // Values of N
      READ( NIN, FMT = * )NIDIM;
      if ( NIDIM < 1 || NIDIM > NIDMAX ) {
         WRITE( NOUT, FMT = 9997 )'N', NIDMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM );
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ) < 0 || IDIM( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX;
            GO TO 230;
         }
      } // 10
      // Values of K
      READ( NIN, FMT = * )NKB;
      if ( NKB < 1 || NKB > NKBMAX ) {
         WRITE( NOUT, FMT = 9997 )'K', NKBMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( KB( I ), I = 1, NKB );
      for (I = 1; I <= NKB; I++) { // 20
         if ( KB( I ) < 0 ) {
            WRITE( NOUT, FMT = 9995 );
            GO TO 230;
         }
      } // 20
      // Values of INCX and INCY
      READ( NIN, FMT = * )NINC;
      if ( NINC < 1 || NINC > NINMAX ) {
         WRITE( NOUT, FMT = 9997 )'INCX AND INCY', NINMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( INC( I ), I = 1, NINC );
      for (I = 1; I <= NINC; I++) { // 30
         if ( INC( I ) == 0 || ( INC( I ) ).abs() > INCMAX ) {
            WRITE( NOUT, FMT = 9994 )INCMAX;
            GO TO 230;
         }
      } // 30
      // Values of ALPHA
      READ( NIN, FMT = * )NALF;
      if ( NALF < 1 || NALF > NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF );
      // Values of BETA
      READ( NIN, FMT = * )NBET;
      if ( NBET < 1 || NBET > NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET );

      // Report values of parameters.

      WRITE( NOUT, FMT = 9993 );
      WRITE( NOUT, FMT = 9992 )( IDIM( I ), I = 1, NIDIM );
      WRITE( NOUT, FMT = 9991 )( KB( I ), I = 1, NKB );
      WRITE( NOUT, FMT = 9990 )( INC( I ), I = 1, NINC );
      WRITE( NOUT, FMT = 9989 )( ALF( I ), I = 1, NALF );
      WRITE( NOUT, FMT = 9988 )( BET( I ), I = 1, NBET );
      if ( !TSTERR ) {
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = 9980 );
      }
      WRITE( NOUT, FMT = * );
      WRITE( NOUT, FMT = 9999 )THRESH;
      WRITE( NOUT, FMT = * );

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      for (I = 1; I <= NSUBS; I++) { // 40
         LTEST[I] = false;
      } // 40
   50 READ( NIN, FMT = 9984, END = 80 )SNAMET, LTESTT;
      for (I = 1; I <= NSUBS; I++) { // 60
         if( SNAMET == SNAMES( I ) ) GO TO 70;
      } // 60
      WRITE( NOUT, FMT = 9986 )SNAMET;
      STOP;
   70 LTEST( I ) = LTESTT;
      GO TO 50;

      } // 80
      CLOSE ( NIN );

      // Compute EPS (the machine precision).

      EPS = EPSILON(RZERO);
      WRITE( NOUT, FMT = 9998 )EPS;

      // Check the reliability of ZMVCH using exact data.

      N = min( 32, NMAX );
      for (J = 1; J <= N; J++) { // 120
         for (I = 1; I <= N; I++) { // 110
            A[I][J] = max( I - J + 1, 0 );
         } // 110
         X[J] = J;
         Y[J] = ZERO;
      } // 120
      for (J = 1; J <= N; J++) { // 130
         YY[J] = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
      } // 130
      // YY holds the exact result. On exit from ZMVCH YT holds
      // the result computed by ZMVCH.
      TRANS = 'N';
      zmvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LZE( YY, YT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR;
         STOP;
      }
      TRANS = 'T';
      zmvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LZE( YY, YT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR;
         STOP;
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 210
         WRITE( NOUT, FMT = * );
         if ( !LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM );
         } else {
            srnamc.SRNAMT = SNAMES( ISNUM );
            // Test error exits.
            if ( TSTERR ) {
               zchke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * );
            }
            // Test computations.
            infoc.INFOT = 0;
            infoc.OK = true;
            FATAL = false;
            GO TO ( 140, 140, 150, 150, 150, 160, 160, 160, 160, 160, 160, 170, 170, 180, 180, 190, 190 )ISNUM;
            // Test ZGEMV, 01, and ZGBMV, 02.
  140       CALL ZCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G );
            GO TO 200;
            // Test ZHEMV, 03, ZHBMV, 04, and ZHPMV, 05.
  150       CALL ZCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G );
            GO TO 200;
            // Test ZTRMV, 06, ZTBMV, 07, ZTPMV, 08,
            // ZTRSV, 09, ZTBSV, 10, and ZTPSV, 11.
  160       CALL ZCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z );
            GO TO 200;
            // Test ZGERC, 12, ZGERU, 13.
  170       CALL ZCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z );
            GO TO 200;
            // Test ZHER, 14, and ZHPR, 15.
  180       CALL ZCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z );
            GO TO 200;
            // Test ZHER2, 16, and ZHPR2, 17.
  190       CALL ZCHK6( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z );

  200       IF( FATAL && SFATAL ) GO TO 220;
         }
      } // 210
      WRITE( NOUT, FMT = 9982 );
      GO TO 240;

      } // 220
      WRITE( NOUT, FMT = 9981 );
      GO TO 240;

      } // 230
      WRITE( NOUT, FMT = 9987 );

      } // 240
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT );
      STOP;

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LESS THAN${.f8_2}');
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE${( * 10).d9_1}');
 9997 FORMAT( ' NUMBER OF VALUES OF ${} IS LESS THAN 1 OR GREATER THAN ${.i2}');
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${.i2}');
 9995 FORMAT( ' VALUE OF K IS LESS THAN 0' );
 9994 FORMAT( ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ${.i2}');
 9993 FORMAT( ' TESTS OF THE Complex       LEVEL 2 BLAS', //' THE FOLLOWING PARAMETER VALUES WILL BE USED:' );
 9992 FORMAT( '   FOR N              ', 9I6 );
 9991 FORMAT( '   FOR K              ', 7I6 );
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 );
 9989 FORMAT( '   FOR ALPHA          ', 7( '(${.f4_1},${.f4_1})  ', : ) );
 9988 FORMAT( '   FOR BETA           ', 7( '(${.f4_1},${.f4_1})  ', : ) );
 9987 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******' );
 9986 FORMAT( ' SUBPROGRAM NAME ${.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******' );
 9985 FORMAT( ' ERROR IN ZMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n ZMVCH WAS CALLED WITH TRANS = ${.a1} AND RETURNED SAME = ${.l1} AND ERR = ${.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******' );
 9984 FORMAT( A6, L2 );
 9983 FORMAT(' ${.a6} WAS NOT TESTED' );
 9982 FORMAT('\n END OF TESTS' );
 9981 FORMAT('\n ******* FATAL ERROR - TESTS ABANDONED *******' );
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' );
      }
      void zchk1(final int SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, final int FATAL, final int NIDIM, final int IDIM, final int NKB, final int KB, final int NALF, final int ALF, final int NBET, final int BET, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int Y, final int YY, final int YS, final int YT, final int G,) {

// Tests ZGEMV and ZGBMV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, BETA, BLS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA, LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK, NL, NS;
      bool               BANDED, FULL, NULL, RESET, SAME, TRAN;
      String             TRANS, TRANSS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZGBMV, ZGEMV, ZMAKE, ZMVCH, ZREGR1
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Data statements ..
      const ICH = 'NTC';
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'E';
      BANDED = SNAME( 3: 3 ) == 'B';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 11;
      } else if ( BANDED ) {
         NARGS = 13;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN );
         ND = N/2 + 1;

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = max( N - ND, 0 );
            IF( IM == 2 ) M = min( N + ND, NMAX );

            if ( BANDED ) {
               NK = NKB;
            } else {
               NK = 1;
            }
            for (IKU = 1; IKU <= NK; IKU++) { // 100
               if ( BANDED ) {
                  KU = KB( IKU );
                  KL = max( KU - 1, 0 );
               } else {
                  KU = N - 1;
                  KL = M - 1;
               }
               // Set LDA to 1 more than minimum value if room.
               if ( BANDED ) {
                  LDA = KL + KU + 1;
               } else {
                  LDA = M;
               }
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 100;
               LAA = LDA*N;
               NULL = N <= 0 || M <= 0;

               // Generate the matrix A.

               TRANSL = ZERO;
               zmake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL );

               for (IC = 1; IC <= 3; IC++) { // 90
                  TRANS = ICH( IC: IC );
                  TRAN = TRANS == 'T' || TRANS == 'C';

                  if ( TRAN ) {
                     ML = N;
                     NL = M;
                  } else {
                     ML = M;
                     NL = N;
                  }

                  for (IX = 1; IX <= NINC; IX++) { // 80
                     INCX = INC( IX );
                     LX = ( INCX ).abs()*NL;

                     // Generate the vector X.

                     TRANSL = HALF;
                     zmake('GE', ' ', ' ', 1, NL, X, 1, XX, ( INCX ).abs(), 0, NL - 1, RESET, TRANSL );
                     if ( NL > 1 ) {
                        X[NL/2] = ZERO;
                        XX[1 + ( INCX ).abs()*( NL/2 - 1 )] = ZERO;
                     }

                     for (IY = 1; IY <= NINC; IY++) { // 70
                        INCY = INC( IY );
                        LY = ( INCY ).abs()*ML;

                        for (IA = 1; IA <= NALF; IA++) { // 60
                           ALPHA = ALF( IA );

                           for (IB = 1; IB <= NBET; IB++) { // 50
                              BETA = BET( IB );

                              // Generate the vector Y.

                              TRANSL = ZERO;
                              zmake('GE', ' ', ' ', 1, ML, Y, 1, YY, ( INCY ).abs(), 0, ML - 1, RESET, TRANSL );

                              NC = NC + 1;

                              // Save every datum before calling the
                              // subroutine.

                              TRANSS = TRANS;
                              MS = M;
                              NS = N;
                              KLS = KL;
                              KUS = KU;
                              ALS = ALPHA;
                              for (I = 1; I <= LAA; I++) { // 10
                                 AS[I] = AA( I );
                              } // 10
                              LDAS = LDA;
                              for (I = 1; I <= LX; I++) { // 20
                                 XS[I] = XX( I );
                              } // 20
                              INCXS = INCX;
                              BLS = BETA;
                              for (I = 1; I <= LY; I++) { // 30
                                 YS[I] = YY( I );
                              } // 30
                              INCYS = INCY;

                              // Call the subroutine.

                              if ( FULL ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 zgemv(TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              } else if ( BANDED ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 zgbmv(TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              }

                              // Check if error-exit was taken incorrectly.

                              if ( !infoc.OK ) {
                                 WRITE( NOUT, FMT = 9993 );
                                 FATAL = true;
                                 GO TO 130;
                              }

                              // See what data changed inside subroutines.

                              ISAME[1] = TRANS == TRANSS;
                              ISAME[2] = MS == M;
                              ISAME[3] = NS == N;
                              if ( FULL ) {
                                 ISAME[4] = ALS == ALPHA;
                                 ISAME[5] = LZE( AS, AA, LAA );
                                 ISAME[6] = LDAS == LDA;
                                 ISAME[7] = LZE( XS, XX, LX );
                                 ISAME[8] = INCXS == INCX;
                                 ISAME[9] = BLS == BETA;
                                 if ( NULL ) {
                                    ISAME[10] = LZE( YS, YY, LY );
                                 } else {
                                    ISAME[10] = LZERES( 'GE', ' ', 1, ML, YS, YY, ( INCY ).abs() );
                                 }
                                 ISAME[11] = INCYS == INCY;
                              } else if ( BANDED ) {
                                 ISAME[4] = KLS == KL;
                                 ISAME[5] = KUS == KU;
                                 ISAME[6] = ALS == ALPHA;
                                 ISAME[7] = LZE( AS, AA, LAA );
                                 ISAME[8] = LDAS == LDA;
                                 ISAME[9] = LZE( XS, XX, LX );
                                 ISAME[10] = INCXS == INCX;
                                 ISAME[11] = BLS == BETA;
                                 if ( NULL ) {
                                    ISAME[12] = LZE( YS, YY, LY );
                                 } else {
                                    ISAME[12] = LZERES( 'GE', ' ', 1, ML, YS, YY, ( INCY ).abs() );
                                 }
                                 ISAME[13] = INCYS == INCY;
                              }

                              // If data was incorrectly changed, report
                              // and return.

                              SAME = true;
                              for (I = 1; I <= NARGS; I++) { // 40
                                 SAME = SAME && ISAME( I );
                                 if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                              } // 40
                              if ( !SAME ) {
                                 FATAL = true;
                                 GO TO 130;
                              }

                              if ( !NULL ) {

                                 // Check the result.

                                 zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                                 ERRMAX = max( ERRMAX, ERR );
                                 // If got really bad answer, report and
                                 // return.
                                 if (FATAL) GO TO 130;
                              } else {
                                 // Avoid repeating tests with M <= 0 or
                                 // N <= 0.
                                 GO TO 110;
                              }

                           } // 50

                        } // 60

                     } // 70

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Regression test to verify preservation of y when m zero, n nonzero.

      zregr1(TRANS, M, N, LY, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY, YS );
      if ( FULL ) {
         if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
         if (REWI) REWIND NTRA;
         zgemv(TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
      } else if ( BANDED ) {
         if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
         if (REWI) REWIND NTRA;
         zgbmv(TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
      }
      NC = NC + 1;
      if ( !LZE( YS, YY, LY ) ) {
         WRITE( NOUT, FMT = 9998 )NARGS - 1;
         FATAL = true;
         GO TO 130;
      }

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 140;

      } // 130
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
      }

      } // 140
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${i3(4, ',')}', '(${.f4_1},${.f4_1}), A,${.i3}, X,${.i2},(${.f4_1},${.f4_1}), Y,${.i2}) .' );
 9994 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${i3(2, ',')}', '(${.f4_1},${.f4_1}), A,${.i3}, X,${.i2},(${.f4_1},${.f4_1}), Y,${.i2})         .' );
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchk2(final int SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, final int FATAL, final int NIDIM, final int IDIM, final int NKB, final int KB, final int NALF, final int ALF, final int NBET, final int BET, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int Y, final int YY, final int YS, final int YT, final int G,) {

// Tests ZHEMV, ZHBMV and ZHPMV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, BETA, BLS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IK, IN, INCX, INCXS, INCY, INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHBMV, ZHEMV, ZHPMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'E';
      BANDED = SNAME( 3: 3 ) == 'B';
      PACKED = SNAME( 3: 3 ) == 'P';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 10;
      } else if ( BANDED ) {
         NARGS = 11;
      } else if ( PACKED ) {
         NARGS = 9;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN );

         if ( BANDED ) {
            NK = NKB;
         } else {
            NK = 1;
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK );
            } else {
               K = N - 1;
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1;
            } else {
               LDA = N;
            }
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2;
            } else {
               LAA = LDA*N;
            }
            NULL = N <= 0;

            for (IC = 1; IC <= 2; IC++) { // 90
               UPLO = ICH( IC: IC );

               // Generate the matrix A.

               TRANSL = ZERO;
               zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

               for (IX = 1; IX <= NINC; IX++) { // 80
                  INCX = INC( IX );
                  LX = ( INCX ).abs()*N;

                  // Generate the vector X.

                  TRANSL = HALF;
                  zmake('GE', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     X[N/2] = ZERO;
                     XX[1 + ( INCX ).abs()*( N/2 - 1 )] = ZERO;
                  }

                  for (IY = 1; IY <= NINC; IY++) { // 70
                     INCY = INC( IY );
                     LY = ( INCY ).abs()*N;

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA );

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB );

                           // Generate the vector Y.

                           TRANSL = ZERO;
                           zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );

                           NC = NC + 1;

                           // Save every datum before calling the
                           // subroutine.

                           UPLOS = UPLO;
                           NS = N;
                           KS = K;
                           ALS = ALPHA;
                           for (I = 1; I <= LAA; I++) { // 10
                              AS[I] = AA( I );
                           } // 10
                           LDAS = LDA;
                           for (I = 1; I <= LX; I++) { // 20
                              XS[I] = XX( I );
                           } // 20
                           INCXS = INCX;
                           BLS = BETA;
                           for (I = 1; I <= LY; I++) { // 30
                              YS[I] = YY( I );
                           } // 30
                           INCYS = INCY;

                           // Call the subroutine.

                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              zhemv(UPLO, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              zhbmv(UPLO, N, K, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              zhpmv(UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( !infoc.OK ) {
                              WRITE( NOUT, FMT = 9992 );
                              FATAL = true;
                              GO TO 120;
                           }

                           // See what data changed inside subroutines.

                           ISAME[1] = UPLO == UPLOS;
                           ISAME[2] = NS == N;
                           if ( FULL ) {
                              ISAME[3] = ALS == ALPHA;
                              ISAME[4] = LZE( AS, AA, LAA );
                              ISAME[5] = LDAS == LDA;
                              ISAME[6] = LZE( XS, XX, LX );
                              ISAME[7] = INCXS == INCX;
                              ISAME[8] = BLS == BETA;
                              if ( NULL ) {
                                 ISAME[9] = LZE( YS, YY, LY );
                              } else {
                                 ISAME[9] = LZERES( 'GE', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME[10] = INCYS == INCY;
                           } else if ( BANDED ) {
                              ISAME[3] = KS == K;
                              ISAME[4] = ALS == ALPHA;
                              ISAME[5] = LZE( AS, AA, LAA );
                              ISAME[6] = LDAS == LDA;
                              ISAME[7] = LZE( XS, XX, LX );
                              ISAME[8] = INCXS == INCX;
                              ISAME[9] = BLS == BETA;
                              if ( NULL ) {
                                 ISAME[10] = LZE( YS, YY, LY );
                              } else {
                                 ISAME[10] = LZERES( 'GE', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME[11] = INCYS == INCY;
                           } else if ( PACKED ) {
                              ISAME[3] = ALS == ALPHA;
                              ISAME[4] = LZE( AS, AA, LAA );
                              ISAME[5] = LZE( XS, XX, LX );
                              ISAME[6] = INCXS == INCX;
                              ISAME[7] = BLS == BETA;
                              if ( NULL ) {
                                 ISAME[8] = LZE( YS, YY, LY );
                              } else {
                                 ISAME[8] = LZERES( 'GE', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME[9] = INCYS == INCY;
                           }

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 40
                              SAME = SAME && ISAME( I );
                              if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                           } // 40
                           if ( !SAME ) {
                              FATAL = true;
                              GO TO 120;
                           }

                           if ( !NULL ) {

                              // Check the result.

                              zmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                              ERRMAX = max( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 120;
                           } else {
                              // Avoid repeating tests with N <= 0
                              GO TO 110;
                           }

                        } // 50

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY;
      }

      } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},(${.f4_1},${.f4_1}), AP, X,${.i2},(${.f4_1},${.f4_1}), Y,${.i2})                .' );
 9994 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${i3(2, ',')}', '(${.f4_1},${.f4_1}), A,${.i3}, X,${.i2},(${.f4_1},${.f4_1}), Y,${.i2})         .' );
 9993 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},(${.f4_1},${.f4_1}), A,${.i3}, X,${.i2},(${.f4_1},${.f4_1}), Y,${.i2})             .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchk3(final int SNAME, final int EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NKB, final int KB, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int XT, final int G, final int Z,) {

// Tests ZTRMV, ZTBMV, ZTPMV, ZTRSV, ZTBSV and ZTPSV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XT( NMAX ), XX( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         TRANSL;
      double             ERR, ERRMAX;
      int                I, ICD, ICT, ICU, IK, IN, INCX, INCXS, IX, K, KS, LAA, LDA, LDAS, LX, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             DIAG, DIAGS, TRANS, TRANSS, UPLO, UPLOS;
      String             ICHD, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZMAKE, ZMVCH, ZTBMV, ZTBSV, ZTPMV, ZTPSV, ZTRMV, ZTRSV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Data statements ..
      const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN';
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'R';
      BANDED = SNAME( 3: 3 ) == 'B';
      PACKED = SNAME( 3: 3 ) == 'P';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 8;
      } else if ( BANDED ) {
         NARGS = 9;
      } else if ( PACKED ) {
         NARGS = 7;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;
      // Set up zero vector for ZMVCH.
      for (I = 1; I <= NMAX; I++) { // 10
         Z[I] = ZERO;
      } // 10

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN );

         if ( BANDED ) {
            NK = NKB;
         } else {
            NK = 1;
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK );
            } else {
               K = N - 1;
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1;
            } else {
               LDA = N;
            }
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2;
            } else {
               LAA = LDA*N;
            }
            NULL = N <= 0;

            for (ICU = 1; ICU <= 2; ICU++) { // 90
               UPLO = ICHU( ICU: ICU );

               for (ICT = 1; ICT <= 3; ICT++) { // 80
                  TRANS = ICHT( ICT: ICT );

                  for (ICD = 1; ICD <= 2; ICD++) { // 70
                     DIAG = ICHD( ICD: ICD );

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 2: 3 ), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

                     for (IX = 1; IX <= NINC; IX++) { // 60
                        INCX = INC( IX );
                        LX = ( INCX ).abs()*N;

                        // Generate the vector X.

                        TRANSL = HALF;
                        zmake('GE', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
                        if ( N > 1 ) {
                           X[N/2] = ZERO;
                           XX[1 + ( INCX ).abs()*( N/2 - 1 )] = ZERO;
                        }

                        NC = NC + 1;

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO;
                        TRANSS = TRANS;
                        DIAGS = DIAG;
                        NS = N;
                        KS = K;
                        for (I = 1; I <= LAA; I++) { // 20
                           AS[I] = AA( I );
                        } // 20
                        LDAS = LDA;
                        for (I = 1; I <= LX; I++) { // 30
                           XS[I] = XX( I );
                        } // 30
                        INCXS = INCX;

                        // Call the subroutine.

                        if ( SNAME( 4: 5 ) == 'MV' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              ztrmv(UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              ztbmv(UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              ztpmv(UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        } else if ( SNAME( 4: 5 ) == 'SV' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              ztrsv(UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              ztbsv(UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              ztpsv(UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !infoc.OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 120;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = UPLO == UPLOS;
                        ISAME[2] = TRANS == TRANSS;
                        ISAME[3] = DIAG == DIAGS;
                        ISAME[4] = NS == N;
                        if ( FULL ) {
                           ISAME[5] = LZE( AS, AA, LAA );
                           ISAME[6] = LDAS == LDA;
                           if ( NULL ) {
                              ISAME[7] = LZE( XS, XX, LX );
                           } else {
                              ISAME[7] = LZERES( 'GE', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME[8] = INCXS == INCX;
                        } else if ( BANDED ) {
                           ISAME[5] = KS == K;
                           ISAME[6] = LZE( AS, AA, LAA );
                           ISAME[7] = LDAS == LDA;
                           if ( NULL ) {
                              ISAME[8] = LZE( XS, XX, LX );
                           } else {
                              ISAME[8] = LZERES( 'GE', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME[9] = INCXS == INCX;
                        } else if ( PACKED ) {
                           ISAME[5] = LZE( AS, AA, LAA );
                           if ( NULL ) {
                              ISAME[6] = LZE( XS, XX, LX );
                           } else {
                              ISAME[6] = LZERES( 'GE', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME[7] = INCXS == INCX;
                        }

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                        } // 40
                        if ( !SAME ) {
                           FATAL = true;
                           GO TO 120;
                        }

                        if ( !NULL ) {
                           if ( SNAME( 4: 5 ) == 'MV' ) {

                              // Check the result.

                              zmvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, true );
                           } else if ( SNAME( 4: 5 ) == 'SV' ) {

                              // Compute approximation to original vector.

                              for (I = 1; I <= N; I++) { // 50
                                 Z[I] = XX( 1 + ( I - 1 )* ( INCX ).abs() )                                  XX( 1 + ( I - 1 )*( INCX ).abs() ) = X( I );
                              } // 50
                              zmvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, false );
                           }
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 120;
                        } else {
                           // Avoid repeating tests with N <= 0.
                           GO TO 110;
                        }

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX;
      }

      } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(', 3( '''${.a1}'',' ), I3, ', AP, X,${.i2})                                      .' );
 9994 FORMAT(' ${.i6}: ${.a6}(', 3( '''${.a1}'',' ), 2( I3, ',' ), ' A,${.i3}, X,${.i2})                               .' );
 9993 FORMAT(' ${.i6}: ${.a6}(', 3( '''${.a1}'',' ), I3, ', A,${.i3}, X,${.i2})                                   .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchk4(final int SNAME, EPS, THRESH, NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int Y, final int YY, final int YS, final int YT, final int G, final int Z,) {

// Tests ZGERC and ZGERU.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS, NC, ND, NS;
      bool               CONJ, NULL, RESET, SAME;
      // .. Local Arrays ..
      Complex         W( 1 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZGERC, ZGERU, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Executable Statements ..
      CONJ = SNAME( 5: 5 ) == 'C';
      // Define the number of arguments.
      NARGS = 9;

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN );
         ND = N/2 + 1;

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = max( N - ND, 0 );
            IF( IM == 2 ) M = min( N + ND, NMAX );

            // Set LDA to 1 more than minimum value if room.
            LDA = M;
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 110;
            LAA = LDA*N;
            NULL = N <= 0 || M <= 0;

            for (IX = 1; IX <= NINC; IX++) { // 100
               INCX = INC( IX );
               LX = ( INCX ).abs()*M;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('GE', ' ', ' ', 1, M, X, 1, XX, ( INCX ).abs(), 0, M - 1, RESET, TRANSL );
               if ( M > 1 ) {
                  X[M/2] = ZERO;
                  XX[1 + ( INCX ).abs()*( M/2 - 1 )] = ZERO;
               }

               for (IY = 1; IY <= NINC; IY++) { // 90
                  INCY = INC( IY );
                  LY = ( INCY ).abs()*N;

                  // Generate the vector Y.

                  TRANSL = ZERO;
                  zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     Y[N/2] = ZERO;
                     YY[1 + ( INCY ).abs()*( N/2 - 1 )] = ZERO;
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 80
                     ALPHA = ALF( IA );

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, M - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1;

                     // Save every datum before calling the subroutine.

                     MS = M;
                     NS = N;
                     ALS = ALPHA;
                     for (I = 1; I <= LAA; I++) { // 10
                        AS[I] = AA( I );
                     } // 10
                     LDAS = LDA;
                     for (I = 1; I <= LX; I++) { // 20
                        XS[I] = XX( I );
                     } // 20
                     INCXS = INCX;
                     for (I = 1; I <= LY; I++) { // 30
                        YS[I] = YY( I );
                     } // 30
                     INCYS = INCY;

                     // Call the subroutine.

                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA;
                     if ( CONJ ) {
                        if (REWI) REWIND NTRA;
                        zgerc(M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     } else {
                        if (REWI) REWIND NTRA;
                        zgeru(M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( !infoc.OK ) {
                        WRITE( NOUT, FMT = 9993 );
                        FATAL = true;
                        GO TO 140;
                     }

                     // See what data changed inside subroutine.

                     ISAME[1] = MS == M;
                     ISAME[2] = NS == N;
                     ISAME[3] = ALS == ALPHA;
                     ISAME[4] = LZE( XS, XX, LX );
                     ISAME[5] = INCXS == INCX;
                     ISAME[6] = LZE( YS, YY, LY );
                     ISAME[7] = INCYS == INCY;
                     if ( NULL ) {
                        ISAME[8] = LZE( AS, AA, LAA );
                     } else {
                        ISAME[8] = LZERES( 'GE', ' ', M, N, AS, AA, LDA );
                     }
                     ISAME[9] = LDAS == LDA;

                     // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME && ISAME( I );
                        if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                     } // 40
                     if ( !SAME ) {
                        FATAL = true;
                        GO TO 140;
                     }

                     if ( !NULL ) {

                        // Check the result column by column.

                        if ( INCX > 0 ) {
                           for (I = 1; I <= M; I++) { // 50
                              Z[I] = X( I );
                           } // 50
                        } else {
                           for (I = 1; I <= M; I++) { // 60
                              Z[I] = X( M - I + 1 );
                           } // 60
                        }
                        for (J = 1; J <= N; J++) { // 70
                           if ( INCY > 0 ) {
                              W[1] = Y( J );
                           } else {
                              W[1] = Y( N - J + 1 );
                           }
                           if (CONJ) W( 1 ) = DCONJG( W( 1 ) );
                           zmvch('N', M, 1, ALPHA, Z, NMAX, W, 1, ONE, A( 1, J ), 1, YT, G, AA( 1 + ( J - 1 )*LDA ), EPS, ERR, FATAL, NOUT, true );
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 130;
                        } // 70
                     } else {
                        // Avoid repeating tests with M <= 0 or N <= 0.
                        GO TO 110;
                     }

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 150;

      } // 130
      WRITE( NOUT, FMT = 9995 )J;

      } // 140
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA;

      } // 150
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(${i3(2, ',')}', '(${.f4_1},${.f4_1}), X,${.i2}, Y,${.i2}, A,${.i3})                         .' );
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchk5(final int SNAME, EPS, THRESH, NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int Y, final int YY, final int YS, final int YT, final int G, final int Z,) {

// Tests ZHER and ZHPR.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex         ALPHA, TRANSL;
      double             ERR, ERRMAX, RALPHA, RALS;
      int                I, IA, IC, IN, INCX, INCXS, IX, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      Complex         W( 1 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHER, ZHPR, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, MAX
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'E';
      PACKED = SNAME( 3: 3 ) == 'P';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 7;
      } else if ( PACKED ) {
         NARGS = 6;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN );
         // Set LDA to 1 more than minimum value if room.
         LDA = N;
         if (LDA < NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA > NMAX) GO TO 100;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2;
         } else {
            LAA = LDA*N;
         }

         for (IC = 1; IC <= 2; IC++) { // 90
            UPLO = ICH( IC: IC );
            UPPER = UPLO == 'U';

            for (IX = 1; IX <= NINC; IX++) { // 80
               INCX = INC( IX );
               LX = ( INCX ).abs()*N;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('GE', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
               if ( N > 1 ) {
                  X[N/2] = ZERO;
                  XX[1 + ( INCX ).abs()*( N/2 - 1 )] = ZERO;
               }

               for (IA = 1; IA <= NALF; IA++) { // 70
                  RALPHA = (ALF( IA )).toDouble();
                  ALPHA = DCMPLX( RALPHA, RZERO );
                  NULL = N <= 0 || RALPHA == RZERO;

                  // Generate the matrix A.

                  TRANSL = ZERO;
                  zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                  NC = NC + 1;

                  // Save every datum before calling the subroutine.

                  UPLOS = UPLO;
                  NS = N;
                  RALS = RALPHA;
                  for (I = 1; I <= LAA; I++) { // 10
                     AS[I] = AA( I );
                  } // 10
                  LDAS = LDA;
                  for (I = 1; I <= LX; I++) { // 20
                     XS[I] = XX( I );
                  } // 20
                  INCXS = INCX;

                  // Call the subroutine.

                  if ( FULL ) {
                     if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, RALPHA, INCX, LDA;
                     if (REWI) REWIND NTRA;
                     zher(UPLO, N, RALPHA, XX, INCX, AA, LDA );
                  } else if ( PACKED ) {
                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, RALPHA, INCX;
                     if (REWI) REWIND NTRA;
                     zhpr(UPLO, N, RALPHA, XX, INCX, AA );
                  }

                  // Check if error-exit was taken incorrectly.

                  if ( !infoc.OK ) {
                     WRITE( NOUT, FMT = 9992 );
                     FATAL = true;
                     GO TO 120;
                  }

                  // See what data changed inside subroutines.

                  ISAME[1] = UPLO == UPLOS;
                  ISAME[2] = NS == N;
                  ISAME[3] = RALS == RALPHA;
                  ISAME[4] = LZE( XS, XX, LX );
                  ISAME[5] = INCXS == INCX;
                  if ( NULL ) {
                     ISAME[6] = LZE( AS, AA, LAA );
                  } else {
                     ISAME[6] = LZERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA );
                  }
                  if ( !PACKED ) {
                     ISAME[7] = LDAS == LDA;
                  }

                  // If data was incorrectly changed, report and return.

                  SAME = true;
                  for (I = 1; I <= NARGS; I++) { // 30
                     SAME = SAME && ISAME( I );
                     if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                  } // 30
                  if ( !SAME ) {
                     FATAL = true;
                     GO TO 120;
                  }

                  if ( !NULL ) {

                     // Check the result column by column.

                     if ( INCX > 0 ) {
                        for (I = 1; I <= N; I++) { // 40
                           Z[I] = X( I );
                        } // 40
                     } else {
                        for (I = 1; I <= N; I++) { // 50
                           Z[I] = X( N - I + 1 );
                        } // 50
                     }
                     JA = 1;
                     for (J = 1; J <= N; J++) { // 60
                        W[1] = DCONJG( Z( J ) );
                        if ( UPPER ) {
                           JJ = 1;
                           LJ = J;
                        } else {
                           JJ = J;
                           LJ = N - J + 1;
                        }
                        zmvch('N', LJ, 1, ALPHA, Z( JJ ), LJ, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
                        if ( FULL ) {
                           if ( UPPER ) {
                              JA = JA + LDA;
                           } else {
                              JA = JA + LDA + 1;
                           }
                        } else {
                           JA = JA + LJ;
                        }
                        ERRMAX = max( ERRMAX, ERR );
                        // If got really bad answer, report and return.
                        if (FATAL) GO TO 110;
                     } // 60
                  } else {
                     // Avoid repeating tests if N <= 0.
                     if (N <= 0) GO TO 100;
                  }

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 110
      WRITE( NOUT, FMT = 9995 )J;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, RALPHA, INCX, LDA;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, RALPHA, INCX;
      }

      } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},${.f4_1}, X,${.i2}, AP)                                         .' );
 9993 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},${.f4_1}, X,${.i2}, A,${.i3})                                      .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchk6(final int SNAME, EPS, THRESH, NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NINC, final int INC, final int NMAX, final int INCMAX, final int A, final int AA, final int AS, final int X, final int XX, final int XS, final int Y, final int YY, final int YS, final int YT, final int G, final int Z,) {

// Tests ZHER2 and ZHPR2.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX, 2 );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      Complex         W( 2 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHER2, ZHPR2, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'E';
      PACKED = SNAME( 3: 3 ) == 'P';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 9;
      } else if ( PACKED ) {
         NARGS = 8;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 140
         N = IDIM( IN );
         // Set LDA to 1 more than minimum value if room.
         LDA = N;
         if (LDA < NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA > NMAX) GO TO 140;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2;
         } else {
            LAA = LDA*N;
         }

         for (IC = 1; IC <= 2; IC++) { // 130
            UPLO = ICH( IC: IC );
            UPPER = UPLO == 'U';

            for (IX = 1; IX <= NINC; IX++) { // 120
               INCX = INC( IX );
               LX = ( INCX ).abs()*N;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('GE', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
               if ( N > 1 ) {
                  X[N/2] = ZERO;
                  XX[1 + ( INCX ).abs()*( N/2 - 1 )] = ZERO;
               }

               for (IY = 1; IY <= NINC; IY++) { // 110
                  INCY = INC( IY );
                  LY = ( INCY ).abs()*N;

                  // Generate the vector Y.

                  TRANSL = ZERO;
                  zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     Y[N/2] = ZERO;
                     YY[1 + ( INCY ).abs()*( N/2 - 1 )] = ZERO;
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 100
                     ALPHA = ALF( IA );
                     NULL = N <= 0 || ALPHA == ZERO;

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1;

                     // Save every datum before calling the subroutine.

                     UPLOS = UPLO;
                     NS = N;
                     ALS = ALPHA;
                     for (I = 1; I <= LAA; I++) { // 10
                        AS[I] = AA( I );
                     } // 10
                     LDAS = LDA;
                     for (I = 1; I <= LX; I++) { // 20
                        XS[I] = XX( I );
                     } // 20
                     INCXS = INCX;
                     for (I = 1; I <= LY; I++) { // 30
                        YS[I] = YY( I );
                     } // 30
                     INCYS = INCY;

                     // Call the subroutine.

                     if ( FULL ) {
                        if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA;
                        if (REWI) REWIND NTRA;
                        zher2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     } else if ( PACKED ) {
                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY;
                        if (REWI) REWIND NTRA;
                        zhpr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA );
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( !infoc.OK ) {
                        WRITE( NOUT, FMT = 9992 );
                        FATAL = true;
                        GO TO 160;
                     }

                     // See what data changed inside subroutines.

                     ISAME[1] = UPLO == UPLOS;
                     ISAME[2] = NS == N;
                     ISAME[3] = ALS == ALPHA;
                     ISAME[4] = LZE( XS, XX, LX );
                     ISAME[5] = INCXS == INCX;
                     ISAME[6] = LZE( YS, YY, LY );
                     ISAME[7] = INCYS == INCY;
                     if ( NULL ) {
                        ISAME[8] = LZE( AS, AA, LAA );
                     } else {
                        ISAME[8] = LZERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA );
                     }
                     if ( !PACKED ) {
                        ISAME[9] = LDAS == LDA;
                     }

                     // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME && ISAME( I );
                        if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                     } // 40
                     if ( !SAME ) {
                        FATAL = true;
                        GO TO 160;
                     }

                     if ( !NULL ) {

                        // Check the result column by column.

                        if ( INCX > 0 ) {
                           for (I = 1; I <= N; I++) { // 50
                              Z[I][1] = X( I );
                           } // 50
                        } else {
                           for (I = 1; I <= N; I++) { // 60
                              Z[I][1] = X( N - I + 1 );
                           } // 60
                        }
                        if ( INCY > 0 ) {
                           for (I = 1; I <= N; I++) { // 70
                              Z[I][2] = Y( I );
                           } // 70
                        } else {
                           for (I = 1; I <= N; I++) { // 80
                              Z[I][2] = Y( N - I + 1 );
                           } // 80
                        }
                        JA = 1;
                        for (J = 1; J <= N; J++) { // 90
                           W[1] = ALPHA*DCONJG( Z( J, 2 ) );
                           W[2] = DCONJG( ALPHA )*DCONJG( Z( J, 1 ) );
                           if ( UPPER ) {
                              JJ = 1;
                              LJ = J;
                           } else {
                              JJ = J;
                              LJ = N - J + 1;
                           }
                           zmvch('N', LJ, 2, ONE, Z( JJ, 1 ), NMAX, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
                           if ( FULL ) {
                              if ( UPPER ) {
                                 JA = JA + LDA;
                              } else {
                                 JA = JA + LDA + 1;
                              }
                           } else {
                              JA = JA + LJ;
                           }
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 150;
                        } // 90
                     } else {
                        // Avoid repeating tests with N <= 0.
                        if (N <= 0) GO TO 140;
                     }

                  } // 100

               } // 110

            } // 120

         } // 130

      } // 140

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 170;

      } // 150
      WRITE( NOUT, FMT = 9995 )J;

      } // 160
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY;
      }

      } // 170
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},(${.f4_1},${.f4_1}), X,${.i2}, Y,${.i2}, AP)                            .' );
 9993 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',${.i3},(${.f4_1},${.f4_1}), X,${.i2}, Y,${.i2}, A,${.i3})                         .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void zchke(ISNUM, srnamc.SRNAMT, NOUT ) {

// Tests the error exits from the Level 2 Blas.
// Requires a special version of the error-handling routine XERBLA.
// ALPHA, RALPHA, BETA, A, X and Y should not need to be defined.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      int                ISNUM, NOUT;
      String             srnamc.SRNAMT;
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUTC;
      bool               infoc.LERR, infoc.OK;
      // .. Local Scalars ..
      Complex         ALPHA, BETA;
      double             RALPHA;
      // .. Local Arrays ..
      Complex         A( 1, 1 ), X( 1 ), Y( 1 );
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGBMV, ZGEMV, ZGERC, ZGERU, ZHBMV, ZHEMV, ZHER, ZHER2, ZHPMV, ZHPR, ZHPR2, ZTBMV, ZTBSV, ZTPMV, ZTPSV, ZTRMV, ZTRSV
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, infoc.NOUTC, infoc.OK, infoc.LERR
      // .. Executable Statements ..
      // infoc.OK is set to false by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      infoc.OK = true;
      // infoc.LERR is set to true by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      infoc.LERR = false;
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170 )ISNUM;
   10 infoc.INFOT = 1;
      zgemv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgemv('N', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgemv('N', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgemv('N', 2, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgemv('N', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      zgemv('N', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   20 infoc.INFOT = 1;
      zgbmv('/', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgbmv('N', -1, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgbmv('N', 0, -1, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgbmv('N', 0, 0, -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgbmv('N', 2, 0, 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgbmv('N', 0, 0, 1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   30 infoc.INFOT = 1;
      zhemv('/', 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zhemv('U', -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zhemv('U', 2, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zhemv('U', 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zhemv('U', 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   40 infoc.INFOT = 1;
      zhbmv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zhbmv('U', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zhbmv('U', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zhbmv('U', 0, 1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   50 infoc.INFOT = 1;
      zhpmv('/', 0, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zhpmv('U', -1, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zhpmv('U', 0, ALPHA, A, X, 0, BETA, Y, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zhpmv('U', 0, ALPHA, A, X, 1, BETA, Y, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   60 infoc.INFOT = 1;
      ztrmv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrmv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztrmv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrmv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztrmv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztrmv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   70 infoc.INFOT = 1;
      ztbmv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztbmv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztbmv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztbmv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztbmv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztbmv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztbmv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   80 infoc.INFOT = 1;
      ztpmv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpmv('U', '/', 'N', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpmv('U', 'N', '/', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpmv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztpmv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
   90 infoc.INFOT = 1;
      ztrsv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrsv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztrsv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrsv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztrsv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztrsv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  100 infoc.INFOT = 1;
      ztbsv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztbsv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztbsv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztbsv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztbsv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztbsv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztbsv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  110 infoc.INFOT = 1;
      ztpsv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpsv('U', '/', 'N', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpsv('U', 'N', '/', 0, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpsv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztpsv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  120 infoc.INFOT = 1;
      zgerc(-1, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgerc(0, -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgerc(0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgerc(0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgerc(2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  130 infoc.INFOT = 1;
      zgeru(-1, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeru(0, -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgeru(0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgeru(0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgeru(2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  140 infoc.INFOT = 1;
      zher('/', 0, RALPHA, X, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zher('U', -1, RALPHA, X, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zher('U', 0, RALPHA, X, 0, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zher('U', 2, RALPHA, X, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  150 infoc.INFOT = 1;
      zhpr('/', 0, RALPHA, X, 1, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zhpr('U', -1, RALPHA, X, 1, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zhpr('U', 0, RALPHA, X, 0, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  160 infoc.INFOT = 1;
      zher2('/', 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zher2('U', -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zher2('U', 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zher2('U', 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zher2('U', 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      GO TO 180;
  170 infoc.INFOT = 1;
      zhpr2('/', 0, ALPHA, X, 1, Y, 1, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zhpr2('U', -1, ALPHA, X, 1, Y, 1, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zhpr2('U', 0, ALPHA, X, 0, Y, 1, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zhpr2('U', 0, ALPHA, X, 1, Y, 0, A );
      chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

  180 IF( infoc.OK )THEN;
         WRITE( NOUT, FMT = 9999 )srnamc.SRNAMT;
      } else {
         WRITE( NOUT, FMT = 9998 )srnamc.SRNAMT;
      }
      return;

 9999 FORMAT( ' ${.a6} PASSED THE TESTS OF ERROR-EXITS' );
 9998 FORMAT( ' ******* ${.a6} FAILED THE TESTS OF ERROR-EXITS *******' );
      }
      void zmake(final int TYPE, final int UPLO, final int DIAG, final int M, final int N, final int A, final int NMAX, final int AA, final int LDA, final int KL, final int KU, final int RESET, final int TRANSL,) {

// Generates values for an M by N matrix A within the bandwidth
// defined by KL and KU.
// Stores the values in the array AA in the data structure required
// by the routine, with unwanted elements set to rogue value.

// TYPE is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e10, 1.0e10 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      double             RROGUE;
      const              RROGUE = -1.0e10 ;
      // .. Scalar Arguments ..
      Complex         TRANSL;
      int                KL, KU, LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex         A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      //- Complex         ZBEG;
      // EXTERNAL ZBEG
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DCONJG, MAX, MIN
      // .. Executable Statements ..
      GEN = TYPE( 1: 1 ) == 'G';
      SYM = TYPE( 1: 1 ) == 'H';
      TRI = TYPE( 1: 1 ) == 'T';
      UPPER = ( SYM || TRI ) && UPLO == 'U';
      LOWER = ( SYM || TRI ) && UPLO == 'L';
      UNIT = TRI && DIAG == 'U';

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN || ( UPPER && I <= J ) || ( LOWER && I >= J ) ) {
               if( ( I <= J && J - I <= KU ) || ( I >= J && I - J <= KL ) ) {
                  A[I][J] = ZBEG( RESET ) + TRANSL;
               } else {
                  A[I][J] = ZERO;
               }
               if ( I != J ) {
                  if ( SYM ) {
                     A[J][I] = DCONJG( A( I, J ) );
                  } else if ( TRI ) {
                     A[J][I] = ZERO;
                  }
               }
            }
         } // 10
         if (SYM) A( J, J ) = DCMPLX( (A( J, J )).toDouble(), RZERO );
         if[TRI ) A( J][J] = A( J, J ) + ONE;
         IF[UNIT ) A( J][J] = ONE;
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA[I + ( J - 1 )*LDA] = A( I, J );
            } // 30
            for (I = M + 1; I <= LDA; I++) { // 40
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 40
         } // 50
      } else if ( TYPE == 'GB' ) {
         for (J = 1; J <= N; J++) { // 90
            for (I1 = 1; I1 <= KU + 1 - J; I1++) { // 60
               AA[I1 + ( J - 1 )*LDA] = ROGUE;
            } // 60
            for (I2 = I1; I2 <= min( KL + KU + 1, KU + 1 + M - J ); I2++) { // 70
               AA[I2 + ( J - 1 )*LDA] = A( I2 + J - KU - 1, J );
            } // 70
            for (I3 = I2; I3 <= LDA; I3++) { // 80
               AA[I3 + ( J - 1 )*LDA] = ROGUE;
            } // 80
         } // 90
      } else if ( TYPE == 'HE' || TYPE == 'TR' ) {
         for (J = 1; J <= N; J++) { // 130
            if ( UPPER ) {
               IBEG = 1;
               if ( UNIT ) {
                  IEND = J - 1;
               } else {
                  IEND = J;
               }
            } else {
               if ( UNIT ) {
                  IBEG = J + 1;
               } else {
                  IBEG = J;
               }
               IEND = N;
            }
            for (I = 1; I <= IBEG - 1; I++) { // 100
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 100
            for (I = IBEG; I <= IEND; I++) { // 110
               AA[I + ( J - 1 )*LDA] = A( I, J );
            } // 110
            for (I = IEND + 1; I <= LDA; I++) { // 120
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 120
            if ( SYM ) {
               JJ = J + ( J - 1 )*LDA;
               AA[JJ] = DCMPLX( (AA( JJ )).toDouble(), RROGUE );
            }
         } // 130
      } else if ( TYPE == 'HB' || TYPE == 'TB' ) {
         for (J = 1; J <= N; J++) { // 170
            if ( UPPER ) {
               KK = KL + 1;
               IBEG = max( 1, KL + 2 - J );
               if ( UNIT ) {
                  IEND = KL;
               } else {
                  IEND = KL + 1;
               }
            } else {
               KK = 1;
               if ( UNIT ) {
                  IBEG = 2;
               } else {
                  IBEG = 1;
               }
               IEND = min( KL + 1, 1 + M - J );
            }
            for (I = 1; I <= IBEG - 1; I++) { // 140
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 140
            for (I = IBEG; I <= IEND; I++) { // 150
               AA[I + ( J - 1 )*LDA] = A( I + J - KK, J );
            } // 150
            for (I = IEND + 1; I <= LDA; I++) { // 160
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 160
            if ( SYM ) {
               JJ = KK + ( J - 1 )*LDA;
               AA[JJ] = DCMPLX( (AA( JJ )).toDouble(), RROGUE );
            }
         } // 170
      } else if ( TYPE == 'HP' || TYPE == 'TP' ) {
         IOFF = 0;
         for (J = 1; J <= N; J++) { // 190
            if ( UPPER ) {
               IBEG = 1;
               IEND = J;
            } else {
               IBEG = J;
               IEND = N;
            }
            for (I = IBEG; I <= IEND; I++) { // 180
               IOFF = IOFF + 1;
               AA[IOFF] = A( I, J );
               if ( I == J ) {
                  if (UNIT) AA( IOFF ) = ROGUE;
                  IF[SYM ) AA( IOFF] = DCMPLX( (AA( IOFF )).toDouble(), RROGUE );
               }
            } // 180
         } // 190
      }
      return;
      }
      void zmvch(final int TRANS, final int M, final int N, final int ALPHA, final int A, final int NMAX, final int X, final int INCX, final int BETA, final int Y, final int INCY, final int YT, final int G, final int YY, final int EPS, final int ERR, final int FATAL, final int NOUT, final int MV,) {

// Checks the results of the computational tests.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0, RONE = 1.0 ;
      // .. Scalar Arguments ..
      Complex         ALPHA, BETA;
      double             EPS, ERR;
      int                INCX, INCY, M, N, NMAX, NOUT;
      bool               FATAL, MV;
      String             TRANS;
      // .. Array Arguments ..
      Complex         A( NMAX, * ), X( * ), Y( * ), YT( * ), YY( * );
      double             G( * );
      // .. Local Scalars ..
      Complex         C;
      double             ERRI;
      int                I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
      bool               CTRAN, TRAN;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, SQRT
      // .. Statement Functions ..
      double             ABS1;
      // .. Statement Function definitions ..
      ABS1[C] = ( C.toDouble() ).abs() + ( DIMAG( C ) ).abs();
      // .. Executable Statements ..
      TRAN = TRANS == 'T';
      CTRAN = TRANS == 'C';
      if ( TRAN || CTRAN ) {
         ML = N;
         NL = M;
      } else {
         ML = M;
         NL = N;
      }
      if ( INCX < 0 ) {
         KX = NL;
         INCXL = -1;
      } else {
         KX = 1;
         INCXL = 1;
      }
      if ( INCY < 0 ) {
         KY = ML;
         INCYL = -1;
      } else {
         KY = 1;
         INCYL = 1;
      }

      // Compute expected result in YT using data in A, X and Y.
      // Compute gauges in G.

      IY = KY;
      for (I = 1; I <= ML; I++) { // 40
         YT[IY] = ZERO;
         G[IY] = RZERO;
         JX = KX;
         if ( TRAN ) {
            for (J = 1; J <= NL; J++) { // 10
               YT[IY] = YT( IY ) + A( J, I )*X( JX );
               G[IY] = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 10
         } else if ( CTRAN ) {
            for (J = 1; J <= NL; J++) { // 20
               YT[IY] = YT( IY ) + DCONJG( A( J, I ) )*X( JX );
               G[IY] = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 20
         } else {
            for (J = 1; J <= NL; J++) { // 30
               YT[IY] = YT( IY ) + A( I, J )*X( JX );
               G[IY] = G( IY ) + ABS1( A( I, J ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 30
         }
         YT[IY] = ALPHA*YT( IY ) + BETA*Y( IY );
         G[IY] = ABS1( ALPHA )*G( IY ) + ABS1( BETA )*ABS1( Y( IY ) );
         IY = IY + INCYL;
      } // 40

      // Compute the error ratio for this result.

      ERR = ZERO;
      for (I = 1; I <= ML; I++) { // 50
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*( INCY ).abs() ) )/EPS;
         if( G( I ) != RZERO ) ERRI = ERRI/G( I );
         ERR = max( ERR, ERRI );
         if( ERR*sqrt( EPS ) >= RONE ) GO TO 60;
      } // 50
      // If the loop completes, all results are at least half accurate.
      GO TO 80;

      // Report fatal error.

   60 FATAL = true;
      WRITE( NOUT, FMT = 9999 );
      for (I = 1; I <= ML; I++) { // 70
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, YT( I ), YY( 1 + ( I - 1 )*( INCY ).abs() );
         } else {
            WRITE( NOUT, FMT = 9998 )I, YY( 1 + ( I - 1 )*( INCY ).abs() ), YT( I );
         }
      } // 70

      } // 80
      return;

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT' );
 9998 FORMAT( 1X, I7, 2( '  (${.g15_6},${.g15_6})' ) );
      }
      bool lze(final int RI, final int RJ, final int LR,) {

// Tests if two arrays are identical.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      int                LR;
      // .. Array Arguments ..
      Complex         RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         if( RI( I ) != RJ( I ) ) GO TO 20;
      } // 10
      LZE = true;
      GO TO 30;
      } // 20
      LZE = false;
   30 return;
      }
      bool lzeres(final int TYPE, final int UPLO, final int M, final int N, final int AA, final int AS, final int LDA,) {

// Tests if selected elements in two arrays are equal.

// TYPE is 'GE', 'HE' or 'HP'.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex         AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO == 'U';
      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 10
         } // 20
      } else if ( TYPE == 'HE' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1;
               IEND = J;
            } else {
               IBEG = J;
               IEND = N;
            }
            for (I = 1; I <= IBEG - 1; I++) { // 30
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 30
            for (I = IEND + 1; I <= LDA; I++) { // 40
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 40
         } // 50
      }

      LZERES = true;
      GO TO 80;
      } // 70
      LZERES = false;
   80 return;
      }
      Complex FUNCTION ZBEG( RESET );

// Generates complex numbers as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ;
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891;
         MJ = 457;
         I = 7;
         J = 7;
         IC = 0;
         RESET = false;
      }

      // The sequence of values of I or J is bounded between 1 and 999.
      // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I or J = 4 or 8, the period will be 25.
      // If initial I or J = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I or J
      // in 6.

      IC = IC + 1;
   10 I = I*MI;
      J = J*MJ;
      I = I - 1000*( I/1000 );
      J = J - 1000*( J/1000 );
      if ( IC >= 5 ) {
         IC = 0;
         GO TO 10;
      }
      ZBEG = DCMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 );
      return;
      }
      double ddiff(final int X, final int Y,) {

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      double             X, Y;
      // .. Executable Statements ..
      DDIFF = X - Y;
      return;
      }
      void chkxer(srnamc.SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK ) {

// Tests whether XERBLA has detected an error when it should.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      int                infoc.INFOT, NOUT;
      bool               infoc.LERR, infoc.OK;
      String             srnamc.SRNAMT;
      // .. Executable Statements ..
      if ( !infoc.LERR ) {
         WRITE( NOUT, FMT = 9999 )infoc.INFOT, srnamc.SRNAMT;
         infoc.OK = false;
      }
      infoc.LERR = false;
      return;

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${.i2} NOT DETECTED BY ${.a6} *****' );
      }
      void zregr1(final int TRANS, final int M, final int N, final int LY, final int KL, final int KU, final int ALPHA, final Matrix<double> A_, final int LDA, final int X, final int INCX, final int BETA, final int Y, final int INCY, final int YS,) {
  final A = A_.dim();

// Input initialization for regression test.
      String             TRANS;
      int                LY, M, N, KL, KU, LDA, INCX, INCY;
      Complex         ALPHA, BETA;
      // .. Array Arguments ..
      Complex         A(LDA,*), X(*), Y(*), YS(*);
      // .. Local Scalars ..
      int                I;
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // .. Executable Statements ..
      TRANS = 'T';
      M = 0;
      N = 5;
      KL = 0;
      KU = 0;
      ALPHA = DCMPLX( 1.0 );
      LDA = max( 1, M );
      INCX = 1;
      BETA = DCMPLX( -0.7, -0.8 );
      INCY = 1;
      LY = ( INCY ).abs()*N;
      for (I = 1; I <= LY; I++) { // 10
         Y[I] = DCMPLX( 42.0, I.toDouble() );
         YS[I] = Y( I );
      } // 10
      return;
      }
      void xerbla(final int SRNAME, final Box<int> INFO,) {

// This is a special version of XERBLA to be used only as part of
// the test program for testing error exits from the Level 2 BLAS
// routines.

// XERBLA  is an error handler for the Level 2 BLAS routines.

// It is called by the Level 2 BLAS routines if an input parameter is
// invalid.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.
      int                INFO;
      String             SRNAME;
      // .. Scalars in Common ..
      int                infoc.INFOT, NOUT;
      bool               infoc.LERR, infoc.OK;
      String             srnamc.SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON /SRNAMC/srnamc.SRNAMT
      // .. Executable Statements ..
      infoc.LERR = true;
      if ( INFO != infoc.INFOT ) {
         if ( infoc.INFOT != 0 ) {
            WRITE( NOUT, FMT = 9999 )INFO, infoc.INFOT;
         } else {
            WRITE( NOUT, FMT = 9997 )INFO;
         }
         infoc.OK = false;
      }
      if ( SRNAME != srnamc.SRNAMT ) {
         WRITE( NOUT, FMT = 9998 )SRNAME, srnamc.SRNAMT;
         infoc.OK = false;
      }
      return;

 9999 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ${.i6} INSTEAD OF ${.i2} *******' );
 9998 FORMAT( ' ******* XERBLA WAS CALLED WITH SRNAME = ${.a6} INSTEAD OF ${.a6} *******' );
 9997 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ${.i6} *******' );
      }
