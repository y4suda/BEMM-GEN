# CEMM-gen の実装

## TODO
* ディレクトリ構造の決定
* 処理の流れを考える
* 必要なライブラリを探す
  * 引数の処理
  * setuptools
  * ドキュメント化
* チュートリアル作成

## 基本的な処理の流れ

## ディレクトリ構造

## ライブラリ選定

### 引数の処理
標準ライブラリの argparse を使用

サブコマンドを分けてパーサーを用意
https://qiita.com/kzkadc/items/e4fc7bc9c003de1eb6d0
https://qiita.com/oohira/items/308bbd33a77200a35a3d


### インストール方法
https://packaging.python.org/ja/latest/guides/distributing-packages-using-setuptools/
https://qiita.com/propella/items/5cd89caee6379920d889

### ドキュメント化
Sphinxがいい？
https://zenn.dev/yamagishihrd/articles/8b5279a07000a4

docstring のスタイルはGoogle style と Numpy styleがある→とりあえず google style?