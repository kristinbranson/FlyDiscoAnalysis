% SET   AVIFILE オブジェクトのプロパティの設定
%
% OBJ = SET(OBJ,'PropertyName',VALUE) は、AVIFILE オブジェクト OBJ の
% プロパティ 'PropertyName' に値 VALUE を設定します。 
%
% OBJ = SET(OBJ,'PropertyName',VALUE,'PropertyName',VALUE,..) は、単一の
% ステートメントで、AVIFILE オブジェクト OBJ の複数のプロパティ値を設定します。
%
% 注意: この関数は、SUBSASGN 用の補助関数で、ユーザの使用を目的とした
%       ものではありません。AVIFILE オブジェクトの適切な値を設定するには
%       構造体記法を使ってください。例えばつぎのようにします。:
%
%       obj.Fps = value;


%   Copyright 1984-2006 The MathWorks, Inc.
