%ADDFRAME  ビデオフレームを AVI ファイルに追加
%
%   AVIOBJ = ADDFRAME(AVIOBJ,FRAME) は、FRAME 内のデータを AVIFILE で
%   作成された AVIOBJ に追加します。FRAME はインデックス付きイメージ (M×N)、
%   または、double か uint8 の精度のトゥルーカラーイメージ (M×N×3) の
%   いずれかになります。FRAME が AVI ファイルの 1 番目のフレームに追加されない
%   場合は、前のフレームの次元との整合性を保たなければなりません。
%
%   AVIOBJ = ADDFRAME(AVIOBJ,FRAME1,FRAME2,FRAME3,...) は、複数のフレームを 
%   AVI ファイルに追加します。
%
%   AVIOBJ = ADDFRAME(AVIOBJ,MOV) は、MATLAB ムービー MOV に含まれるフレーム
%   を AVI ファイルに追加します。インデックス付きイメージとしてフレームを
%   格納する MATLAB ムービーは、カラーマップがあらかじめ設定されていない場合、
%   AVI ファイル用のカラーマップとして 1 番目のフレームのカラーマップを使用します。
%
%   AVIOBJ = ADDFRAME(AVIOBJ,H) は、Figure、または、軸のハンドル H から
%   フレームをキャプチャし、このフレームを AVI ファイルに追加します。
%   フレームは、AVI ファイルに追加される前に、画面外の配列に描画されます。
%   アニメーション内のグラフィックスが XOR グラフィックスを使用している場合、
%   このシンタックスを使用する必要はありません。
%
%   アニメーションが XOR グラフィックスを使用している場合、MATLAB ムービーの 
%   1 つのフレームにグラフィックスをキャプチャする代わりに GETFRAME を使用し、
%   下記の例のように、[AVIOBJ] = ADDFRAME(AVIOBJ,MOV) を使用します。
%   GETFRAME は、画面上のイメージのスナップショットを撮ります。
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
%   参考 AVIFILE, AVIFILE/CLOSE, MOVIE2AVI.


%   Copyright 1984-2008 The MathWorks, Inc.
