function time = tocf(mode) 
% TOCF 時間の表示を成形したストップウオッチタイマの読み込み
%
%	TOCF(M)
%	TIC と TOCF 関数は、経過時間を測定するために一緒に機能します。
%	TOC 自身では、TIC を使用してからの経過時間を整数値Mで指定された形式で表示します。
%	t = TOCF; は、経過日数、時間、分、秒の要素を持つ数値配列として t に経過時間を保存します。
%
%	表示形式については下記の表を参考にしてください。
%	引数として表示形式を与えなかった場合は，1～4の範囲で0の要素がないものだけが表示される形式が選択されます。
%
%	---------------------------------------------------------------
%	| 整数値 |      表示形式
%	---------------------------------------------------------------
%	|    1   |   経過時間は11.123456秒です 
%	|    2   |   経過時間は12分11.123456秒です 
%	|    3   |   経過時間は1時間12分11.123456秒です 
%	|    4   |   経過時間は1日1時間12分11.123456秒です 
%	|    5   |   Elapsed time is 12.123456 seconds. 
%	|    6   |   Elapsed time is 12 min 11.123456 sec. 
%	|    7   |   Elapsed time is 1 hr 12 min 11.123456 sec. 
%	|    8   |   Elapsed time is 1 day 1 hr 12 min 11.123456 sec. 
%	---------------------------------------------------------------
%
%	例
%	ストップウォッチの経過時間を分と秒で表示します。
%	tocf(4)
%
%	参考 toc, tic, cputime.

% --
%	Title : TOCF()
%	Author : Sach1o : http://sach1o.blog80.fc2.com/
%	Created : 2007/12/12
% //--

% 入力チェックと引数なしの場合の処理
if nargin==1 & [ mode<1 | mode>8 | ~isequal(mode, floor(mode))] 
	error('第1引数は1～4の正の整数値である必要があります。');
end;

% 日、時、分、秒に分ける
ntime=toc;
nday = floor(ntime/86400);
nhour = floor(mod(ntime,86400)/3600);
nmin = floor(mod(ntime,3600)/60);
nsec = mod(ntime,60);

% 以下出力処理
if nargout==1 
	time = [nday, nhour, nmin, nsec];
elseif nargout==0 
	if nargin==0 
		if nday>0 
			mode=4;
		elseif nhour>0 
			mode=3;
		elseif nmin>0 
			mode=2;
		else 
			mode=1;
		end;
	end;

	switch mode 
		case 1 
			strtime = sprintf( '%0.6f%s', ntime, '秒');
		case 2 
			nmin = nmin + 60*nhour + 1440*nday;
			strtime = sprintf( '%02d%s%04.6f%s', nmin, '分', nsec, '秒');
		case 3 
			nhour = nhour + 24*nday;
			strtime = sprintf( '%d%s%02d%s%04.6f%s', nhour, '時間',nmin, '分', nsec, '秒');
		case 4 
			strtime = sprintf( '%d%s%d%s%02d%s%04.6f%s', nday, '日', nhour, '時間',nmin, '分', nsec, '秒');
		case 5 
			strtime = sprintf( '%0.6f%s', ntime, ' seconds');
		case 6 
			nmin = nmin + 60*nhour + 1440*nday;
			strtime = sprintf( '%02d%s%04.6f%s', nmin, ' min ', nsec, ' sec');
		case 7 
			nhour = nhour + 24*nday;
			strtime = sprintf( '%d%s%02d%s%04.6f%s', nhour, ' hr ',nmin, ' min ', nsec, ' sec');
		case 8	 
			strtime = sprintf( '%d%s%d%s%02d%s%04.6f%s', nday, ' day ', nhour, ' hr ',nmin, ' min ', nsec, ' sec');
	end;

	if mode<5 
		disp(sprintf( '%s%s%s', '経過時間は', strtime, 'です'));
	else 
		disp(sprintf( '%s%s%s', 'Elapsed time is ', strtime, '.'));
	end;
end;
