# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'SingleStepProblem'
        db.create_table(u'app_singlestepproblem', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reactantSmiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('productSmiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('correctAnswer', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'app', ['SingleStepProblem'])


    def backwards(self, orm):
        # Deleting model 'SingleStepProblem'
        db.delete_table(u'app_singlestepproblem')


    models = {
        u'app.singlestepproblem': {
            'Meta': {'object_name': 'SingleStepProblem'},
            'correctAnswer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'productSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactantSmiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['app']