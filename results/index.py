from flask import Flask, render_template
import pandas
import os

app = Flask(__name__)

def get_thresholds():
    return [("0","0.0"),("1","1.0"),("2","2.0"),("3","3.0"),("4","4.0"),
           ("5","5.0"),("6","6.0"),("7","7.0"),("8","8.0"),("9","9.0")]

def get_metrics():
    return {"pearson":"Pearson","spearman":"Spearman"}

def get_directions():
    return {"posneg":"Positive and Negative","pos":"Positive Only"}

def get_strategies():
    return {"cca":"Complete Case Analysis","svi":"Single Value Imputation"}

def get_accuracy(strategy,metric,direction,threshold):
   df = pandas.read_csv("static/data/supp_data2_ml_accuracy.csv")
   df = df[df.direction==direction]
   df = df[df.thresh==int(threshold)]
   df = df[df.strategy=="%s.%s" %(strategy,metric)]
   return "%.3f (%.3f,%.3f)" %(df.accuracy_mean.tolist()[0],df.down.tolist()[0],df.up.tolist()[0])

@app.route('/about')
def show_about():
    return render_template('about.html')

@app.route('/')
def show_best_analysis():

    metric = "pearson"
    threshold = "1"
    strategy = "cca"
    direction = "posneg"
    datafile = "%s_%s_%s_%s.tsv" %(strategy,metric,direction,threshold)

    # We will show accuracy values
    acc = get_accuracy(strategy,metric,direction,threshold)

    return render_template('index.html',datafile=datafile,
                            strategies=get_strategies(),
                            thresholds=get_thresholds(),
                            metrics=get_metrics(),
                            directions=get_directions(),
                            threshold=threshold,
                            direction=direction,
                            metric=metric,
                            strategy=strategy,
                            acc=acc)

@app.route('/<strategy>/<metric>/<direction>/<threshold>')
def show_specific_analysis(strategy,direction,metric,threshold):

    datafile = "%s_%s_%s_%s.tsv" %(strategy,metric,direction,threshold)

    # We will show accuracy values
    acc = get_accuracy(strategy,metric,direction,threshold)

    return render_template('index.html',datafile=datafile,
                            strategies=get_strategies(),
                            thresholds=get_thresholds(),
                            metrics=get_metrics(),
                            directions=get_directions(),
                            threshold=threshold,
                            direction=direction,
                            strategy=strategy,
                            metric=metric)

if __name__ == '__main__':
    app.debug = True
    app.run()
