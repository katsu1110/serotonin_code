# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:07:15 2019

@author: katsuhisa
"""

# spike prediction from LFP by a simple 1D CNN


# load mat file ===========================================
import scipy.io as sio

mat_contents = sio.loadmat('filepath')

# generate datasets =====================
stlfp1 = np.random.random((1000, 101))
stlfp0 = np.random.random((1000, 101))

X = np.vstack((stlfp1, stlfp0))
y = np.vstack((np.ones((1000,), dtype=int), np.zeros((1000,), dtype=int))

# 1D CNN model ===============================================
model_m = Sequential()
model_m.add(Reshape((TIME_PERIODS, num_sensors), input_shape=(input_shape,)))
model_m.add(Conv1D(100, 10, activation='relu', input_shape=(TIME_PERIODS, num_sensors)))
model_m.add(Conv1D(100, 10, activation='relu'))
model_m.add(MaxPooling1D(3))
model_m.add(Conv1D(160, 10, activation='relu'))
model_m.add(Conv1D(160, 10, activation='relu'))
model_m.add(GlobalAveragePooling1D())
model_m.add(Dropout(0.5))
model_m.add(Dense(num_classes, activation='softmax'))
print(model_m.summary())

# cross-validation (leave-one-out) ======================
callbacks_list = [
    keras.callbacks.ModelCheckpoint(
        filepath='best_model.{epoch:02d}-{val_loss:.2f}.h5',
        monitor='val_loss', save_best_only=True),
    keras.callbacks.EarlyStopping(monitor='acc', patience=1)
]

model_m.compile(loss='binary_crossentropy',
                optimizer='rmsprop', metrics=['accuracy'])

BATCH_SIZE = 400
EPOCHS = 50

history = model_m.fit(x_train,
                      y_train,
                      batch_size=BATCH_SIZE,
                      epochs=EPOCHS,
                      callbacks=callbacks_list,
                      validation_split=0.2,
                      verbose=1)