%AVIFILE  新しい AVI ファイルの作成
%
%   AVIOBJ = AVIFILE(FILENAME) は、デフォルトのパラメータ値を持つ AVIFILE 
%   オブジェクト AVIOBJ を作成します。FILENAME が拡張子を含んでいない場合、
%   '.avi' が使われます。AVIFILE で開かれたファイルを閉じるには、AVIFILE/CLOSE を
%   使用します。開いているすべての AVI ファイルを閉じるには、"clear mex" を使用
%   してください。
%
%   AVIOBJ = AVIFILE(FILENAME,'PropertyName',VALUE,'PropertyName',VALUE,...) 
%   は、指定したプロパティ値を持つ AVIFILE オブジェクトを返します。
%
%   AVIFILE パラメータ
%
%   FPS         - AVI ムービー用の秒毎のフレーム。このパラメータは、ADDFRAME 
%   の使用前に設定されていなければなりません。デフォルトは、15 fps です。
%
%   COMPRESSION - 圧縮に用いる方法を指定する文字列です。UNIX では、この値は 
%   'None' でなければなりません。Windows 用の利用可能なパラメータは、
%   'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE' または 'None' です。
%
%   ユーザの圧縮法を使用する場合、値は符号化のドキュメントで指定された 4 つの
%   文字コードを使用します。指定されたユーザの圧縮法が見つからない場合、
%   ADDFRAME を呼び出している間にエラーになります。このパラメータは、ADDFRAME の
%   使用前に設定されていなければなりません。
%   デフォルトは、Windows では 'Indeo5'、UNIX では 'None' です。
%   注意: Indeo5 は、Windows のいくつかのバージョンで利用できない可能性があります。
%
%   QUALITY      - 0 から 100 の間の数。このパラメータは圧縮されていない
%   ムービーには影響しません。このパラメータは、ADDFRAME の使用前に設定されて
%   いなければなりません。高い値は、高いビデオ画質でより大きなファイルサイズになり、
%   低い値は、低いビデオ画質でより小さなファイルサイズになります。
%   デフォルトは 75 です。
%
%   KEYFRAME     - 一時的な圧縮をサポートする圧縮法に対して、このパラメータ
%   は、単位時間 (秒) あたりのキーフレーム数です。このパラメータは、ADDFRAME の
%   使用前に設定されていなければなりません。デフォルトは毎秒 2 フレームです。
%
%   COLORMAP     - インデックス付き AVI ムービーに対して使用されるカラーマップ
%   を定義する M 行 3 列の行列です。M は 256 より小さい数でなければなりません 
%   (Indeo 圧縮は 236 です)。このパラメータは、MATLAB のムービーシンタックスと
%   して ADDFRAME を使う以外では、ADDFRAME が呼び出される前に設定されていなければ
%   なりません。デフォルトのカラーマップはありません。
%
%   VIDEONAME    - ビデオストリーム用の記述名。このパラメータは、64 文字より
%   小さく、ADDFRAME を使用する前に設定されていなければなりません。デフォルトは
%   ファイル名です。
%
%
%   AVIFILE プロパティは、MATLAB 構造体シンタックスを用いて設定できます。
%   たとえば、以下のシンタックスを使って Quality プロパティに 100 を設定します。
%
%      aviobj = avifile(filename);
%      aviobj.Quality = 100;
%
%
%   例:
%
%      t = linspace(0,2.5*pi,40);
%      fact = 10*sin(t);
%      fig=figure;
%      aviobj = avifile('example.avi')
%      [x,y,z] = peaks;
%      for k=1:length(fact)
%          h = surf(x,y,fact(k)*z);
%          axis([-3 3 -3 3 -80 80])
%          axis off
%          caxis([-90 90])
%          F = getframe(fig);
%          aviobj = addframe(aviobj,F);
%      end
%      close(fig)
%      aviobj = close(aviobj);
%
%
%   参考 AVIFILE/ADDFRAME, AVIFILE/CLOSE, MOVIE2AVI.


%   Copyright 1984-2008 The MathWorks, Inc.
