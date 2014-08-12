# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'PredictProductsProblem'
        db.create_table(u'app_predictproductsproblem', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reactantSmiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('reagents', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('correctAnswer', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'app', ['PredictProductsProblem'])

        # Adding model 'Synthesis'
        db.create_table(u'app_synthesis', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reactantSmiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('productSmiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('correctAnswer', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'app', ['Synthesis'])


    def backwards(self, orm):
        # Deleting model 'PredictProductsProblem'
        db.delete_table(u'app_predictproductsproblem')

        # Deleting model 'Synthesis'
        db.delete_table(u'app_synthesis')


    models = {
        u'app.predictproductsproblem': {
            'Meta': {'object_name': 'PredictProductsProblem'},
            'correctAnswer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reactantSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reagents': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.singlestepproblem': {
            'Meta': {'object_name': 'SingleStepProblem'},
            'correctAnswer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'productSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactantSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.synthesis': {
            'Meta': {'object_name': 'Synthesis'},
            'correctAnswer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'productSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactantSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['app']