# チュートリアル

## インストール方法
```sh
conda create cemm-gen-env
conda activate cemm-gen-env
conda install -c conda-forge openbabel psi4
pip install -e ../CEMM-GEN
```

## 基本的な使い方
```sh:available_sub-commands
# 円筒形のモデルを作成
cemm-gen cylinder

# 球形のモデルを作成
cemm-gen sphere

# 壁面の残基を作成
cemm-gen makeparam

# 使用可能な残基一覧を表示
cemm-gen listparam
```

## 円筒形のモデル作成

## 球形のモデル作成



## タンパク質構造の入力または生成
CEMM-GEN では生成したモデル内にタンパク質を配置できます。タンパク質の指定方法は２種類です。

1. pdbファイルを指定する `--proteinpdb`

    任意のタンパク質構造を含む構造ファイルを指定できます

2. アミノ酸配列と二次構造を指定する `--proteinseq` および `--proteinSS`

    Ambertools の tleap プログラムを使用して構造を自動生成します。アミノ酸配列は１文字表記で入力してください。
    二次構造は直鎖（コイル） `C` もしくはヘリックス `H` のみが使用できます。アミノ酸配列と長さが一致しない場合はエラーが返されます。

いずれの場合も、タンパク質の長径、短径は自動で計算されモデルの大きさに反映されます。`cylinder` の `--length` を指定した場合タンパク質がはみ出す場合は警告が出されますが、構造生成は許容されます。最終的な構造をよく確認してください。`--radius` を指定した場合、タンパク質の大きさがはみ出すことは許容されないため、エラーを返し終了します。タンパク質のサイズに合わせたモデルを作るには `--padding-radius` を活用してください。

## パラメーターの作成

現在利用可能な残基一覧は `listparam` サブコマンドで取得できます。

```sh:quick_example
cemm-gen listparam
```

独自の残基を作る場合は `makeparam` サブコマンドを使用してください。

```sh:quick_example
cemm-gen makeparam --smiles CCC --resname MTY --description "Methyl group"
```

### 最低限必要な入力
- SMILES表記の化合物構造 `--smiles`
  - 位置拘束をかけるため、必ず２つの炭素で開始してください。
- 残基名 `--resname`
  - すでに利用されているアミノ酸残基名と被らないように、３文字で名前をつけてください。
- 説明 `--description`
  - 残基を登録する際の説明を入力してください。

### 構造最適化を省略する場合
Rdkit の ETKDGv3 法で構造生成を行った後、`--method-opt` と `--basisSet-opt` で指定された汎関数？で最適化を行います。
これを省略する場合（構造最適化が完了しないなど）、`--singlePoint` オプションを付けて実行してください。

### 電荷をもつ残基の場合
`--netcharge` オプションに電荷を指定してください。
負電荷を持つカルボキシ基などは芳香族結合 `:` をSMILESに明記してください。 例：`CCC(:O):O`

### 実行速度に関する設定
`--num-thread` に使用するスレッド数を指定してください。`--memory-sizeGB` に使用可能なメモリサイズを指定してください。いずれもマシンの全コア、メモリを使用すると他の操作ができなくなる場合があるので、適宜余裕をもって設定してください。